Param(
  [Parameter(Mandatory=$true)] [string]$VolumeName,
  [Parameter(Mandatory=$true)] [int]$CPUs,
  [Parameter(Mandatory=$true)] [int]$MemoryGB,
  [Parameter(Mandatory=$true)] [string]$ExportDir
)

Write-Host "Volumen   : $VolumeName"
Write-Host "CPUs/RAM  : $CPUs / ${MemoryGB}GB"
Write-Host "ExportDir : $ExportDir"
Write-Host "==> QC de VCFs filtrados (Mutect2) …"

# Bash que se ejecutará en Docker
$bash = @'
set -euo pipefail
set -x
shopt -s nullglob
shopt -s extglob

trap 'rc=$?; echo "[ERROR] comando falló con rc=${rc}"; exit ${rc}' ERR

echo "[env] shell=$SHELL"
REF="/work/ref/GRCh38_primary.fa"
echo "[ref] ${REF}"

OUT="/work"
if [[ "${EXPORT:-}" == "1" ]]; then
  OUT="/out"
fi
mkdir -p "${OUT}"

summary="${OUT}/summary.tsv"
errors="${OUT}/errors.tsv"
echo -e "file\tsample\ttotal_variants\tpass_variants\tpass_rate\tpass_snps\tpass_indels\tts\ttv\ttstv" > "${summary}"
echo -e "file\terror" > "${errors}"

bcftools --version

echo "[find] buscando /work/**/*.somatic.filtered.vcf.gz (recursivo)"
mapfile -t VCF_LIST < <(find /work -type f -name "*.somatic.filtered.vcf.gz" | sort)
echo "[find] candidatos: ${#VCF_LIST[@]}"

for vcf in "${VCF_LIST[@]}"; do
  echo "[QC] ${vcf}"
  base="${vcf%.somatic.filtered.vcf.gz}"
  stats="${base}.somatic.filtered.bcfstats"

  # 1)stats
  if ! bcftools stats -s - "${vcf}" > "${stats}" 2> >(tee /dev/stderr); then
    echo -e "${vcf}\tbcftools stats falló" >> "${errors}"
    continue
  fi

  # 2)totales y PASS
  if ! total=$(bcftools view -H "${vcf}" | wc -l | awk '{print $1}'); then
    echo -e "${vcf}\tbcftools view total falló" >> "${errors}"
    continue
  fi
  if ! pass=$(bcftools view -f PASS -H "${vcf}" | wc -l | awk '{print $1}'); then
    echo -e "${vcf}\tbcftools view PASS falló" >> "${errors}"
    continue
  fi
  pass_rate="0"
  if [[ "${total}" -gt 0 ]]; then
    pass_rate=$(awk -v p="${pass}" -v t="${total}" 'BEGIN{printf("%.6f", p/t)}')
  fi

  # 3)SNPs/INDELs PASS
  if ! pass_snps=$(bcftools view -f PASS -H -v snps "${vcf}" | wc -l | awk '{print $1}'); then
    echo -e "${vcf}\tbcftools view PASS snps falló" >> "${errors}"
    continue
  fi
  if ! pass_indels=$(bcftools view -f PASS -H -v indels "${vcf}" | wc -l | awk '{print $1}'); then
    echo -e "${vcf}\tbcftools view PASS indels falló" >> "${errors}"
    continue
  fi

  # 4)Ts/Tv desde bcftools stats (robusto, sin 'read' ni asumir 'all')
  ts="0"; tv="0"; tstv="NA"
  if grep -q "^TSTV" "${stats}"; then
    ts=$(awk '$1=="TSTV"{print $3; exit}' "${stats}")
    tv=$(awk '$1=="TSTV"{print $4; exit}' "${stats}")
    tstv=$(awk '$1=="TSTV"{print $5; exit}' "${stats}")
    ts=${ts:-0}; tv=${tv:-0}; tstv=${tstv:-"NA"}
  fi

  # 5)sample (primero en el VCF)
  sample=$(bcftools query -l "${vcf}" | head -n1)
  sample=${sample:-"NA"}

  echo -e "${vcf}\t${sample}\t${total}\t${pass}\t${pass_rate}\t${pass_snps}\t${pass_indels}\t${ts}\t${tv}\t${tstv}" >> "${summary}"
done

set +x
echo
echo "[DONE] Resumen QC -> ${summary}"
head -n 10 "${summary}" || true
echo "[NOTE] Errores (si los hubo) -> ${errors}"
'@

# Guardar run.sh como UTF-8 (sin BOM) y con LF
$tmpFile = Join-Path $env:TEMP ("run_" + [System.Guid]::NewGuid().ToString() + ".sh")
$bashUnix = $bash -replace "`r`n","`n"
$enc = New-Object System.Text.UTF8Encoding($false)  # sin BOM
[System.IO.File]::WriteAllText($tmpFile, $bashUnix, $enc)

# Preparar Docker
$image = 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'

$dockerArgs = @(
  'run','--rm','-i',
  '--cpus', "$CPUs",
  '--memory', "${MemoryGB}g",
  '-u','0:0',
  '-v', "${VolumeName}:/work",
  '-v', "${tmpFile}:/root/run.sh:ro"
)

if ($ExportDir -and $ExportDir.Trim() -ne '') {
  if (-not (Test-Path $ExportDir)) { New-Item -ItemType Directory -Path $ExportDir | Out-Null }
  $dockerArgs += @('-v', "${ExportDir}:/out", '-e', 'EXPORT=1')
}

$cmdInContainer = @('/bin/bash','/root/run.sh')

try {
  & docker @dockerArgs $image @cmdInContainer
  $exit = $LASTEXITCODE
} finally {
  if (Test-Path $tmpFile) { Remove-Item $tmpFile -Force }
}

if ($exit -ne 0) {
  throw "falló docker run (ExitCode=$exit)"
} else {
  Write-Host "Listo."
}
