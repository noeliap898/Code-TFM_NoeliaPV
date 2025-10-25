Param(
  [Parameter(Mandatory=$true)] [string]$VolumeName,   #'wesdata'
  [int]$CPUs = 8,
  [int]$MemoryGB = 16
)

Write-Host "Volumen   : $VolumeName"
Write-Host "CPUs/RAM  : $CPUs / ${MemoryGB}GB"
Write-Host "==> Marcando duplicados con GATK (stdout a consola)…"

#Preparamos un script bash dentro del contenedor (sin BOM y con \n)
$runBash = @'
set -euo pipefail

echo "[env] bash=$(bash --version | head -n1)"
echo "[setup] apt-get update && install samtools"
export DEBIAN_FRONTEND=noninteractive
echo 'Acquire::ForceIPv4 "true";' > /etc/apt/apt.conf.d/99force-ipv4
apt-get update -qq
apt-get install -y -qq samtools > /dev/null

mkdir -p /work/tmp /work/aln

echo "[find] buscando BAMs: /work/aln/*.sorted.bam"
shopt -s nullglob
BAMS=(/work/aln/*.sorted.bam)

if (( ${#BAMS[@]} == 0 )); then
  echo "[ERROR] no hay BAMs *.sorted.bam en /work/aln"
  exit 3
fi

echo "[info] encontrados: ${#BAMS[@]}"

for BAM in "${BAMS[@]}"; do
  base="$(basename "$BAM")"
  sample="${base%%.sorted.bam}"
  OUT="/work/aln/${sample}.dedup.bam"
  MET="/work/aln/${sample}.dedup.metrics.txt"
  FLG="/work/aln/${sample}.dedup.flagstat.txt"

  if [[ -f "$OUT" ]]; then
    echo "[skip] ya existe $OUT → salto"
    continue
  fi

  echo "[MarkDuplicates] $base → $(basename "$OUT")"
  gatk --java-options "-Xms2g -Xmx6g" MarkDuplicates \
    -I "$BAM" \
    -O "$OUT" \
    -M "$MET" \
    --CREATE_INDEX true \
    --ASSUME_SORTED true \
    --VALIDATION_STRINGENCY LENIENT \
    --TMP_DIR /work/tmp

  echo "[flagstat] $(basename "$OUT")"
  samtools flagstat "$OUT" > "$FLG"

  echo "[ok] $OUT (+ .bai, metrics y flagstat)"
done

echo "[done] MarkDuplicates terminado"
'@

# Escribimos run.sh sin BOM y con saltos LF
$tmp = Join-Path $env:TEMP "run_markdup.sh"
$utf8NoBom = New-Object System.Text.UTF8Encoding($false)
[IO.File]::WriteAllText($tmp, ($runBash -replace "`r`n","`n"), $utf8NoBom)

# Ejecutamos dentro del contenedor GATK
# solo montamos el volumen con los datos; no hace falta exportar nada ahora.
$dockerArgs = @(
  'run','--rm','--name','wes-dedup',
  '--cpus', "$CPUs",
  '--memory', "${MemoryGB}g",
  '-v', "${VolumeName}:/work",
  '-v', "${tmp}:/root/run.sh:ro",
  'broadinstitute/gatk:4.5.0.0',
  'bash','-lc','bash /root/run.sh'
)

# Si existía de antes, lo intentamos quitar por si quedó colgado
try { docker rm -f wes-dedup *>$null } catch {}

docker @dockerArgs
if ($LASTEXITCODE -ne 0) {
  throw "docker run terminó con código $LASTEXITCODE"
}

Write-Host "Listo. Dedup BAMs y .bai/.flagstat generados en: volumen '$VolumeName' → /work/aln/"
