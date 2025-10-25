param(
  [string]$VolumeName = "wesdata",
  [int]$CPUs = 8,
  [int]$MemoryGB = 16,
  [string]$ExportDir = "D:/results/wes"
)

Write-Host "Volumen   : $VolumeName"
Write-Host "CPUs/RAM  : $CPUs / ${MemoryGB}GB"
Write-Host "ExportDir : $ExportDir"
Write-Host "==> Lanzando Mutect2 sobre BAMs deduplicados (stdout a consola)…"

$bash = @'
set -euo pipefail 2>/dev/null || set -euo

echo "[env] bash=$(bash --version | head -n 1) shell=$SHELL"
Ref=/work/ref/GRCh38_primary.fa
echo "[ref] $Ref"

mkdir -p /work/m2

# Solo organoide WES deduplicados (no RNA-seq)
shopt -s nullglob
found=0
for tumor in /work/aln/*organoid*WES*.dedup.bam; do
  [ -e "$tumor" ] || continue
  found=1

  # ID de paciente = antes del primer '~' del basename
  base=$(basename "$tumor")
  pid="${base%%~*}"

  # Buscar normal por patrón: <PID>~*germlineWES.dedup.bam
  matches=(/work/aln/"${pid}"~*germlineWES.dedup.bam)
  if [ "${#matches[@]}" -ne 1 ] || [ ! -f "${matches[0]}" ]; then
    echo "[WARN] No encontré normal único para $tumor (PID=$pid). Candidates=${#matches[@]}"
    for m in "${matches[@]}"; do echo "  - $m"; done
    continue
  fi
  normal_bam="${matches[0]}"

  tumorSM="${base%.dedup.bam}"
  normalSM="$(basename "$normal_bam" .dedup.bam)"

  echo "[pair] tumor=$tumor"
  echo "[pair] normal=$normal_bam"
  echo "[SM] tumorSM=$tumorSM normalSM=$normalSM"

  outvcf="/work/m2/${tumorSM}.somatic.unfiltered.vcf.gz"

  gatk Mutect2 \
    -R "$Ref" \
    -I "$tumor" -tumor "$tumorSM" \
    -I "$normal_bam" -normal "$normalSM" \
    -O "$outvcf"

  # Indexar por si no se creó
  [ -f "$outvcf" ] && [ ! -f "${outvcf}.tbi" ] && gatk IndexFeatureFile -I "$outvcf" || true
done

if [ "$found" -eq 0 ]; then
  echo "[WARN] No se encontraron BAMs deduplicados de organoide WES en /work/aln/*organoid*WES*.dedup.bam"
fi

mkdir -p /out
cp -v /work/m2/*.vcf.gz /out/ 2>/dev/null || true
cp -v /work/m2/*.vcf.gz.tbi /out/ 2>/dev/null || true

echo "Mutect2-dedup finalizado. VCFs en /work/m2 y copiados a /out."
'@

# Normaliza a LF y guarda sin BOM
$bash = $bash -replace "`r`n", "`n"
$tmp = Join-Path $env:TEMP "run_mutect2.sh"
[IO.File]::WriteAllText($tmp, $bash, (New-Object System.Text.UTF8Encoding($false)))

# Ejecutar contenedor
docker run --rm `
  --name wes-mutect2-dedup `
  -v "${VolumeName}:/work" `
  -v "${ExportDir}:/out" `
  -v "${tmp}:/root/run.sh:ro" `
  broadinstitute/gatk:4.5.0.0 `
  bash -lc "/bin/bash /root/run.sh"
