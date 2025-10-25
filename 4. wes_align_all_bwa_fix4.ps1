#requires -Version 5.1
param(
  [Parameter(Mandatory=$true)][string]$VolumeName,   #'wesdata'
  [Parameter(Mandatory=$true)][int]$CPUs,            #8
  [Parameter(Mandatory=$true)][int]$MemoryGB,        #28
  [Parameter(Mandatory=$true)][string]$ExportDir     #'D:\results\wes'
)

Write-Host ("Volumen   : {0}" -f $VolumeName)
Write-Host ("CPUs/RAM  : {0} / {1}GB" -f $CPUs,$MemoryGB)
Write-Host ("ExportDir : {0}" -f $ExportDir)
New-Item -ItemType Directory -Force -Path $ExportDir | Out-Null
try { docker rm -f wes-all | Out-Null } catch {}

#Script bash a ejecutar dentro del contenedor
$bashScript = @'
#!/usr/bin/env bash
set -eu  # sin pipefail para máxima compatibilidad

echo "[env] bash=${BASH_VERSION:-unknown} shell=$SHELL"
echo 'Acquire::ForceIPv4 "true";' >/etc/apt/apt.conf.d/99force-ipv4
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y bwa samtools gzip

Ref="/work/ref/GRCh38_primary.fa"; [ -f "$Ref" ] || Ref="/work/ref/GRCh38.fa"
echo "[ref] $Ref"
[ -f "$Ref" ] || { echo "[ERROR] no existe $Ref"; exit 20; }

# Índices BWA
if [ ! -f "${Ref}.bwt" ] || [ ! -f "${Ref}.sa" ]; then
  echo "[ERROR] faltan índices BWA (${Ref}.bwt / ${Ref}.sa)"; exit 21
fi

mkdir -p /work/tmp /work/aln
shopt -s nullglob

echo "[find] buscando R1 WES en /work/fastq"
mapfile -t R1s < <(ls -1t \
  /work/fastq/*WES*.R1.*val_1.fq.gz \
  /work/fastq/*WES*.R1.*.fq.gz \
  /work/fastq/*WES*.R1.*.fastq.gz 2>/dev/null || true)
echo "[info] nº R1 detectados: ${#R1s[@]}"
[ "${#R1s[@]}" -gt 0 ] || { echo "[ERROR] no encuentro R1 WES"; exit 3; }

for R1 in "${R1s[@]}"; do
  base="$(basename "$R1")"
  S="${base%%.R1.*}"
  R2="${R1/.R1./.R2.}"; [[ "$R1" == *val_1* ]] && R2="${R2/val_1/val_2}"
  if [ ! -f "$R2" ]; then echo "[WARN] Falta R2 para $base → salto"; continue; fi

  BAM="/work/aln/${S}.sorted.bam"
  if [ -f "$BAM" ]; then echo "[skip] ya existe $BAM"; continue; fi

  echo "[align] muestra=$S"
  echo "[pair]  R1=$R1"
  echo "[pair]  R2=$R2"
  RG="@RG\tID:${S}\tSM:${S}\tPL:ILLUMINA\tLB:WES"

  bwa mem -t "${CPUS:-4}" -K 200000000 -R "$RG" "$Ref" "$R1" "$R2" \
    | samtools sort -@ 2 -m 1G -T /work/tmp -o "$BAM" -

  samtools index -@ 2 "$BAM"
  samtools flagstat -@ 2 "$BAM" > "/work/aln/${S}.flagstat.txt"
  echo "[ok] ${BAM} (+ .bai y flagstat)"
done

echo "[done] alineación para todas las parejas WES"
'@

# Crear fichero temporal LF y montarlo
$tempPath = Join-Path $env:TEMP "wes_align_all.run.sh"
$bashLF = $bashScript -replace "`r`n", "`n"
[IO.File]::WriteAllText($tempPath, $bashLF, [Text.Encoding]::UTF8)

# Ejecutar dentro del contenedor usando el script montado
$runArgs = @(
  'run','--rm',
  '--name=wes-all',
  "--cpus=$CPUs",
  "--memory=${MemoryGB}g",
  '-v',("${VolumeName}:/work"),
  '-v',("${tempPath}:/root/run.sh:ro"),
  '--env',("CPUS=$CPUs"),
  'ubuntu:22.04','/bin/bash','-lc','bash /root/run.sh'
)
& docker @runArgs
if ($LASTEXITCODE -ne 0) { throw "docker run (alineación) terminó con código $LASTEXITCODE" }

Write-Host "Alineación OK. Resultados en el volumen: /work/aln/" -ForegroundColor Green

# Exportar a carpeta del host
Write-Host ("Exportando a {0}..." -f $ExportDir) -ForegroundColor Cyan
$copyArgs = @(
  'run','--rm',
  '-v',("${VolumeName}:/work"),
  '-v',("${ExportDir}:/out"),
  'alpine','sh','-lc',
  'set -e; mkdir -p /out; echo "[ls] /work/aln:"; ls -l /work/aln || true; cp -v /work/aln/* /out/ 2>/dev/null || true; echo "Export hecho."'
)
& docker @copyArgs
if ($LASTEXITCODE -ne 0) {
  Write-Warning "Export falló (¿no hay ficheros en /work/aln?). Revisa el [ls] de arriba."
} else {
  Write-Host "Export completado." -ForegroundColor Green
}

# Limpieza del script temporal
try { Remove-Item $tempPath -ErrorAction SilentlyContinue } catch {}
