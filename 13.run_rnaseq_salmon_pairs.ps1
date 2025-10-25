param(
  [Parameter(Mandatory=$true)] [string]$VolumeName,                 #'wesdata'
  [int]$CPUs = 8,
  [string]$ExportDir = "D:/results/rna/salmon",                     # destino en Windows para copiar resultados
  [string]$IndexDirInVol = "/work/salmon_index",                    # carpeta del índice dentro del volumen
  [string]$TranscriptFaInVol = ""                                   # opcional: ruta EXACTA al FASTA de transcritos dentro del volumen. Si vacío, se autodetecta.
)

Write-Host "Volumen   : $VolumeName"
Write-Host "CPUs      : $CPUs"
Write-Host "Index dir : $IndexDirInVol"
if ($ExportDir) { Write-Host "ExportDir : $ExportDir" }

if ($ExportDir -and !(Test-Path $ExportDir)) {
  New-Item -ItemType Directory -Force -Path $ExportDir | Out-Null
}

# Script bash que corre dentro del contenedor de Salmon
$sh = @'
set -euo pipefail

THREADS="${THREADS:-8}"
IDX_DIR="${IDX_DIR:-/work/salmon_index}"
OUT_BASE="/work/salmon"
RNA_DIR="/work/rna"

TRANS_FA="${TRANS_FA:-}"

echo "==[salmon]== threads=$THREADS"
echo "Index dir : $IDX_DIR"
echo "RNA dir   : $RNA_DIR"
echo "Out base  : $OUT_BASE"

#Validaciones básicas
[ -d "$RNA_DIR" ] || { echo "[ERROR] No existe $RNA_DIR en el volumen (¿copiaste los FASTQ?)"; exit 2; }
mkdir -p "$OUT_BASE"

#Localizar FASTA de transcritos (si no lo pasaron)
if [ -z "$TRANS_FA" ]; then
  echo "Buscando FASTA de transcritos en el volumen (Funcotator GENCODE)…"
  # candidatos típicos dentro de tus data sources (hg38):
  # /work/refs/funcotator/funcotator_dataSources*/gencode/hg38/gencode.v*.pc_transcripts.fa
  TRANS_FA="$(ls -1 /work/refs/funcotator/funcotator_dataSources*/gencode/hg38/gencode.v*.pc_transcripts.fa 2>/dev/null | head -n1 || true)"
fi

if [ -z "$TRANS_FA" ] || [ ! -s "$TRANS_FA" ]; then
  echo "[ERROR] No encuentro el FASTA de transcritos. Pásalo con -TranscriptFaInVol o deja los Funcotator data sources montados en /work/refs/funcotator/."; exit 3;
fi

echo "Transcripts FASTA: $TRANS_FA"

#Construir índice si no existe
if [ ! -s "$IDX_DIR/complete_ref_lens.bin" ]; then
  echo "==[index]== Creando índice Salmon en $IDX_DIR"
  mkdir -p "$IDX_DIR"
  salmon index -t "$TRANS_FA" -i "$IDX_DIR" -p "$THREADS"
else
  echo "==[index]== Ya existe un índice válido en $IDX_DIR (skip)."
fi

#Enumerar parejas R1/R2 y cuantificar
echo "==[find]== Buscando parejas (~RNASEQ.R1/R2.FASTQ.gz) en $RNA_DIR"
mapfile -d '' -t R1S < <(find "$RNA_DIR" -type f -iname '*~rnaseq.r1.fastq.gz' -print0 | sort -z)
echo "R1 detectados: ${#R1S[@]}"

if [ ${#R1S[@]} -eq 0 ]; then
  echo "[WARN] No se encontraron R1 con patrón '~RNASEQ.R1.FASTQ.gz' en $RNA_DIR"; exit 0;
fi

for r1 in "${R1S[@]}"; do
  r2="$(echo "$r1" | sed -E 's/~RNASEQ\.R1\.FASTQ\.gz$/~RNASEQ.R2.FASTQ.gz/I')"
  b1="$(basename "$r1")"
  b2="$(basename "$r2")"

  if [ ! -f "$r2" ]; then
    echo "[SKIP] Falta R2 para: $b1  (esperado: $b2)"
    continue
  fi

  sample="$(basename "$r1" | sed -E 's/~RNASEQ\.R1\.FASTQ\.gz$//I')"
  outdir="$OUT_BASE/$sample"
  if [ -s "$outdir/quant.sf" ]; then
    echo "[skip] Ya existe cuantificación: $outdir/quant.sf"
    continue
  fi

  echo "==[quant]== $sample"
  echo "  R1: $b1"
  echo "  R2: $b2"
  mkdir -p "$outdir"

  salmon quant \
    -i "$IDX_DIR" \
    -l A \
    -1 "$r1" -2 "$r2" \
    -p "$THREADS" \
    --validateMappings \
    --gcBias --seqBias \
    -o "$outdir"
done

#Listado final
echo "==[done]== Resultados (primeros directorios):"
ls -lh "$OUT_BASE" | sed -n '1,60p'
'@

# Guardar script bash temporal con LF
$tmp = [System.IO.Path]::GetTempFileName() + ".sh"
[System.IO.File]::WriteAllText($tmp, ($sh -replace "`r`n","`n"), [System.Text.Encoding]::ASCII)

# Preparar docker args
$dockerArgs = @(
  'run','--rm',
  '--cpus', "$CPUs",
  '-v', "${VolumeName}:/work",
  '-v', "${tmp}:/root/run.sh:ro",
  'combinelab/salmon:1.10.1',     # imagen oficial de Salmon
  'bash','/root/run.sh'
)

# Pasar variables de entorno
$dockerArgs = @(
  'run','--rm',
  '--cpus', "$CPUs",
  '-e', "THREADS=$CPUs",
  '-e', "IDX_DIR=$IndexDirInVol",
  '-v', "${VolumeName}:/work",
  '-v', "${tmp}:/root/run.sh:ro"
)

# Si especificamos FASTA dentro del volumen
if ($TranscriptFaInVol) {
  $dockerArgs += @('-e', "TRANS_FA=$TranscriptFaInVol")
}

# Añadir imagen + comando
$dockerArgs += @('combinelab/salmon:1.10.1','bash','/root/run.sh')

# Ejecutar Salmon
$proc = Start-Process -FilePath 'docker' -ArgumentList $dockerArgs -NoNewWindow -PassThru -Wait
$code = $proc.ExitCode
Remove-Item $tmp -ErrorAction SilentlyContinue

if ($code -ne 0) {
  throw "docker run (salmon) falló (ExitCode=$code)"
}

# Copiar resultados a Windows (si se indicó ExportDir)
if ($ExportDir) {
  Write-Host "==[export]== Copiando /work/salmon -> $ExportDir"
  $copyCmd = @"
set -e
mkdir -p /out
# Copia todo el árbol de resultados
cp -a /work/salmon/. /out/
echo "Contenido exportado:"
ls -lh /out | sed -n '1,120p'
"@
  $tmp2 = [System.IO.Path]::GetTempFileName() + ".sh"
  [System.IO.File]::WriteAllText($tmp2, ($copyCmd -replace "`r`n","`n"), [System.Text.Encoding]::ASCII)

  $args2 = @(
    'run','--rm',
    '-v', "${VolumeName}:/work",
    '-v', "${ExportDir}:/out",
    '-v', "${tmp2}:/root/copy.sh:ro",
    'alpine','sh','/root/copy.sh'
  )
  $p2 = Start-Process -FilePath 'docker' -ArgumentList $args2 -NoNewWindow -PassThru -Wait
  $code2 = $p2.ExitCode
  Remove-Item $tmp2 -ErrorAction SilentlyContinue
  if ($code2 -ne 0) { throw "Exportación falló (ExitCode=$code2)" }
}

Write-Host "Listo."
