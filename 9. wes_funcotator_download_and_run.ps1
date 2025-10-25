param(
  [Parameter(Mandatory = $true )] [string] $VolumeName,              #'wesdata'
  [Parameter(Mandatory = $true )] [int]    $CPUs,                    #8
  [Parameter(Mandatory = $true )] [int]    $MemoryGB,                #24
  [Parameter(Mandatory = $true )] [string] $RefVersion,              #'hg38'
  [Parameter(Mandatory = $true )] [string] $ExportDir,               #'D:/results/wes /funcotator'
  [Parameter(Mandatory = $false)] [string] $DataTarPath = "D:/refs/funcotator/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz",
  [Parameter(Mandatory = $false)] [string] $DataURL     = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz"
)

Write-Host ""
Write-Host "Volumen   : $VolumeName"
Write-Host "CPUs/RAM  : $CPUs / ${MemoryGB}GB"
Write-Host "RefVer    : $RefVersion"
Write-Host "ExportDir : $ExportDir"
Write-Host "DataTar   : $DataTarPath"
Write-Host "DataURL   : $DataURL"
Write-Host "==> Descarga (si hace falta) + anotación Funcotator dentro de Docker…"

#Crear carpetas
$null = New-Item -ItemType Directory -Force -Path (Split-Path $DataTarPath) -ErrorAction SilentlyContinue
$null = New-Item -ItemType Directory -Force -Path $ExportDir -ErrorAction SilentlyContinue

#Descargar con curl.exe si no existe el tar
if (!(Test-Path $DataTarPath)) {
  Write-Host "==> Descargando bundle Funcotator (somático hg38) con curl (reanudable)…"
  $curlArgs = @(
    "-C","-","--fail","-L",
    "-o", "$DataTarPath",
    "$DataURL"
  )
  $proc = Start-Process -FilePath "curl.exe" -ArgumentList $curlArgs -NoNewWindow -PassThru -Wait
  if ($proc.ExitCode -ne 0) {
    throw "curl falló (ExitCode=$($proc.ExitCode)) al descargar $DataURL"
  }
} else {
  Write-Host "==> Ya existe el tar: $DataTarPath  (si estaba parcial, puedes reanudar manualmente con:  curl -C - -o <path> <url>)"
}

#Script bash a ejecutar dentro del contenedor (usa variables de entorno en vez de sed)
$bash = @'
set -euo pipefail
echo "[env] shell=${SHELL:-/bin/bash}  refver=${REFVER:-NA}  xmx=${JAVAMEM_GB:-NA}g"

REF="/work/ref/GRCh38_primary.fa"
FADICT="/work/ref/GRCh38_primary.dict"
FAI="/work/ref/GRCh38_primary.fa.fai"
DATADIR="/opt/funcotator"
DST="$DATADIR/funcotator_dataSources"
TAR="/in/funcotator_dataSources.tar.gz"
OUT="/work"
if [[ "${EXPORT:-}" == "1" ]]; then OUT="/out"; fi
mkdir -p "$OUT"

echo "[check] referencia:"
ls -l "$REF" || true
ls -l "$FADICT" || true
ls -l "$FAI" || true
if [[ ! -s "$REF" || ! -s "$FADICT" || ! -s "$FAI" ]]; then
  echo "[ERROR] Falta FASTA y/o índices (.dict/.fai) en /work/ref/. Esperaba GRCh38_primary.fa + .dict + .fai"
  exit 1
fi

echo "[data] preparando data sources…"
if [[ -d "$DST/hg38" ]]; then
  echo "[data] ya existe $DST/hg38 (saltando extracción)"
else
  if [[ ! -s "$TAR" ]]; then
    echo "[ERROR] No encuentro el tar.gz montado en $TAR"
    exit 1
  fi
  mkdir -p "$DATADIR"
  echo "[untar] extrayendo a $DATADIR (puede tardar)…"
  tar -xzf "$TAR" -C "$DATADIR"
  # Normaliza ruta: a veces queda DATADIR/funcotator_dataSources/hg38
  if [[ -d "$DATADIR/funcotator_dataSources/hg38" && ! -d "$DST/hg38" ]]; then
    ln -s "$DATADIR/funcotator_dataSources" "$DST"
  fi
fi

if [[ ! -d "$DST/hg38" ]]; then
  echo "[ERROR] No veo $DST/hg38 tras extraer."
  ls -la "$DATADIR" || true
  exit 1
fi
echo "[ok] data sources en $DST"

echo "[find] buscando /work/**/*.somatic.filtered.vcf.gz"
mapfile -t VCF_LIST < <(find /work -type f -name "*.somatic.filtered.vcf.gz" | sort)
echo "[find] candidatos: ${#VCF_LIST[@]}"

count=0
for vcf in "${VCF_LIST[@]}"; do
  echo "[Funcotator] $vcf"
  base="${vcf%.somatic.filtered.vcf.gz}"
  bn="$(basename "$base")"
  out="$OUT/${bn}.somatic.funcotator.maf"

  gatk --java-options "-Xmx${JAVAMEM_GB}g" Funcotator \
    --variant "$vcf" \
    --reference "$REF" \
    --ref-version "${REFVER}" \
    --data-sources-path "$DST" \
    --output "$out" \
    --output-file-format MAF
  ((count++)) || true
done

echo "[DONE] MAFs generados: $count"
if [[ "${EXPORT:-}" == "1" ]]; then
  echo "[OUT] Escritos en /out"
fi
'@

# Forzar LF para evitar $'\r'
$bash = $bash -replace "`r`n", "`n"
$tmp = New-TemporaryFile
Set-Content -Path $tmp -Value $bash -NoNewline -Encoding Ascii

#Memoria Java: deja ~4 GB al SO
$JavaGB = [Math]::Max(4, $MemoryGB - 4)

#docker args (ojo con ':', usamos ${}!)
$dockerArgs = @(
  'run','--rm',
  '--cpus', "$CPUs",
  '--memory', "${MemoryGB}g",
  '-e', "JAVAMEM_GB=$JavaGB",
  '-e', "REFVER=$RefVersion",
  '-v', "${VolumeName}:/work",
  '-v', "${tmp}:/root/run.sh:ro",
  '-v', "${DataTarPath}:/in/funcotator_dataSources.tar.gz:ro"
)

if ($ExportDir) {
  $dockerArgs += @('-v', "${ExportDir}:/out", '-e', 'EXPORT=1')
}

$dockerArgs += @(
  'broadinstitute/gatk:4.5.0.0',
  '/bin/bash','-lc',
  'bash /root/run.sh'
)

#Ejecutar
$proc = Start-Process -FilePath 'docker' -ArgumentList $dockerArgs -NoNewWindow -PassThru -Wait
$code = $proc.ExitCode

#Limpieza tmp
Remove-Item $tmp -ErrorAction SilentlyContinue

if ($code -ne 0) {
  throw "falló docker run (ExitCode=$code)"
} else {
  Write-Host "Listo."
}
