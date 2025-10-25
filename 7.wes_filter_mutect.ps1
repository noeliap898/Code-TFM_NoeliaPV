param(
  [string]$VolumeName,
  [int]$CPUs = 4,
  [int]$MemoryGB = 8,
  [string]$ExportDir = "",
  # Ruta de la referencia DENTRO del contenedor
  [string]$ReferenceInside = "/work/ref/GRCh38_primary.fa"
)

Write-Host "Volumen   : $VolumeName"
Write-Host "CPUs/RAM  : $CPUs / ${MemoryGB}GB"
Write-Host "ExportDir : $ExportDir"
Write-Host "==> Lanzando FilterMutectCalls sobre VCFs unfiltered…"

# bash que se ejecuta en el contenedor
$bash = @'
#!/bin/bash
# No usamos -e para que fallos “opcionales” no rompan el contenedor
set -uo pipefail
shopt -s nullglob globstar

REF="${REF:?Referencia no definida (usar -e REF=/ruta/a/ref.fa)}"

echo "[env] bash=$(bash --version | head -1) shell=$SHELL"
echo "[ref] $REF"

if [[ ! -s "$REF" ]]; then
  echo "[ERROR] No encuentro la referencia en $REF" >&2
  # error duro: salimos con 2
  exit 2
fi

echo "[find] buscando /work/**/*.unfiltered.vcf.gz (recursivo)"
mapfile -t VCFs < <(printf "%s\n" /work/**/*.unfiltered.vcf.gz)
echo "[find] candidatos: ${#VCFs[@]}"

processed=0
errors=0

for V in "${VCFs[@]}"; do
  echo "[FilterMutectCalls] -V $V"

  base=${V%.unfiltered.vcf.gz}
  out="${base}.filtered.vcf.gz"

  # stats: probamos nombre estándar; si no, probamos variante sin '.vcf'
  stats1="${base}.unfiltered.vcf.gz.stats"
  stats2="${base}.unfiltered.gz.stats"
  stats=""

  if [[ -s "$stats1" ]]; then
    stats="$stats1"
  elif [[ -s "$stats2" ]]; then
    stats="$stats2"
  fi

  if [[ -n "$stats" ]]; then
    gatk FilterMutectCalls -R "$REF" -V "$V" -O "$out" --stats "$stats" || {
      echo "[ERROR] FilterMutectCalls falló para $V" >&2
      ((errors++))
      continue
    }
  else
    echo "[WARN] no encontré stats para $V — ejecuto sin --stats"
    gatk FilterMutectCalls -R "$REF" -V "$V" -O "$out" || {
      echo "[ERROR] FilterMutectCalls falló para $V (sin stats)" >&2
      ((errors++))
      continue
    }
  fi

  # Indexado (mejor esfuerzo)
  if command -v tabix >/dev/null 2>&1; then
    tabix -f -p vcf "$out" || {
      echo "[WARN] tabix falló, intento gatk IndexFeatureFile"
      gatk IndexFeatureFile -I "$out" || echo "[WARN] también falló IndexFeatureFile"
    }
  else
    echo "[INFO] tabix no disponible; uso gatk IndexFeatureFile"
    gatk IndexFeatureFile -I "$out" || echo "[WARN] IndexFeatureFile falló"
  fi

  # Copiar a /out si está montado (mejor esfuerzo)
  if [[ -n "${EXPORT:-}" ]]; then
    cp -f "$out"* /out/ 2>/dev/null || echo "[WARN] no pude copiar $out a /out (¿montado?)"
  fi

  ((processed++))
done

echo "[DONE] VCFs filtrados: $processed (errores: $errors)"
# No queremos propagar errores no críticos al host: exit 0 siempre que no sea error de referencia
exit 0
'@

# guardar bash con LF
$tmp = New-TemporaryFile
$bashLF = $bash -replace "`r`n", "`n"
[System.IO.File]::WriteAllText($tmp.FullName, $bashLF, (New-Object System.Text.UTF8Encoding($false)))
$tmpPath = $tmp.FullName -replace '\\','/'

# args docker
$dockerArgs = @(
  'run','--rm',
  '--cpus', "$CPUs",
  '--memory', "${MemoryGB}g",
  '-e', "REF=$ReferenceInside",
  '-v', "$($VolumeName):/work",
  '-v', "$($tmpPath):/root/run.sh:ro"
)

if ($ExportDir -and $ExportDir.Trim() -ne '') {
  $dockerArgs += @('-v', "$($ExportDir):/out", '-e', 'EXPORT=1')
}

$dockerArgs += @(
  'broadinstitute/gatk:4.5.0.0',
  'bash','/root/run.sh'
)

# ejecutar docker
& docker @dockerArgs
$exit = $LASTEXITCODE
if ($exit -ne 0) {
  throw "falló docker run (ExitCode=$exit)"
}

Write-Host "Listo"
