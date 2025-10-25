param(
  [Parameter(Mandatory=$true )] [string]$VolumeName,            
  [ValidateSet('hg38','hg19')] [string]$RefVersion = 'hg38',
  [int]$CPUs = 8,
  [int]$MemoryGB = 24,
  [string]$ExportDir = ""                                      
)

Write-Host ""
Write-Host "Volumen   : $VolumeName"
Write-Host "RefVer    : $RefVersion"
Write-Host "CPUs/RAM  : $CPUs / ${MemoryGB}GB"
if ($ExportDir) { Write-Host "ExportDir : $ExportDir" } else { Write-Host "ExportDir : (no se copiará al host)" }

# Memoria Java dejando ~4GB al SO del contenedor
$JavaGB = [Math]::Max(4, $MemoryGB - 4)
if ($ExportDir) { New-Item -ItemType Directory -Force -Path $ExportDir | Out-Null }

function Invoke-DockerScript {
  param(
    [string]$Volume,[string]$Image,[string]$ScriptBody,
    [string]$ExtraMount = "", [hashtable]$Env = @{},
    [int]$CPUs_ = $CPUs, [int]$MemGB_ = $MemoryGB
  )
  $bodyLF = $ScriptBody -replace "`r`n","`n"
  $tmpSh  = [System.IO.Path]::GetTempFileName() + ".sh"
  [System.IO.File]::WriteAllText($tmpSh, $bodyLF, [System.Text.Encoding]::ASCII)

  $args = @('run','--rm','--cpus',"$CPUs_",'--memory',"$($MemGB_)g",
            '-v',"$Volume`:/work",
            '-v',"$tmpSh`:/root/run.sh:ro")
  if ($ExtraMount) { $args += @('-v', $ExtraMount) }
  foreach ($k in $Env.Keys) { $args += @('-e', "$k=$($Env[$k])") }
  $args += @($Image,'bash','/root/run.sh')

  $p = Start-Process -FilePath 'docker' -ArgumentList $args -NoNewWindow -PassThru -Wait
  $code = $p.ExitCode
  Remove-Item $tmpSh -ErrorAction SilentlyContinue
  if ($code -ne 0) { throw "docker run falló (ExitCode=$code)" }
}

# Anotación solo *.filtered.vcf.gz con re-entrada inteligente 
$annotSh = @'
set -euo pipefail

IN_DIR=/work/m2
OUT_DIR=/work/funcotator_out

# Autodetectar raíz del bundle Funcotator
DS_BASE=""
for cand in \
  "/work/refs/funcotator/funcotator_dataSources" \
  /work/refs/funcotator/funcotator_dataSources.* \
  /work/refs/funcotator \
  /work/refs \
; do
  if [ -d "${cand}" ]; then
    n=$(find "${cand}" -maxdepth 3 -type d -name "${REFVER}" 2>/dev/null | wc -l | tr -d " ")
    if [ "${n}" -ge 1 ]; then DS_BASE="${cand}"; break; fi
  fi
done
[ -n "${DS_BASE}" ] || { echo "[ERROR] No pude localizar la raíz del bundle Funcotator (${REFVER})."; exit 3; }

# Referencia
REF_FA=/work/ref/GRCh38.fa
[ "${REFVER}" = "hg19" ] && REF_FA=/work/ref/GRCh37.fa

echo "==[check]=="
echo "REF: ${REF_FA}"
ls -l "${REF_FA}" || { echo "[ERROR] Falta ${REF_FA}"; exit 2; }
[ -f "${REF_FA}.fai" ] || { echo "[ERROR] Falta ${REF_FA}.fai"; exit 2; }
DICT="${REF_FA%.fa}.dict"; [ -f "${DICT}" ] || DICT="${REF_FA}.dict"
[ -f "${DICT}" ] || { echo "[ERROR] Falta .dict de la referencia"; exit 2; }

echo "DS_BASE: ${DS_BASE}"
find "${DS_BASE}" -maxdepth 2 -type d -name "${REFVER}" | head -n 5 || true

echo "IN: ${IN_DIR}"
[ -d "${IN_DIR}" ] || { echo "[ERROR] No existe ${IN_DIR} dentro del volumen"; exit 4; }

echo "==[find]== Buscando *.filtered.vcf.gz en ${IN_DIR}"
mapfile -t VCF_LIST < <(find "${IN_DIR}" -type f -name "*.filtered.vcf.gz" | sort)
echo "Candidatos: ${#VCF_LIST[@]}"
[ ${#VCF_LIST[@]} -gt 0 ] || { echo "[ERROR] No se encontraron VCFs .filtered"; exit 5; }

mkdir -p "${OUT_DIR}"

# Indexar entradas si falta .tbi
echo "==[index]== Creando .tbi si falta"
for V in "${VCF_LIST[@]}"; do
  [ -f "${V}.tbi" ] || tabix -p vcf "${V}" || true
done

echo "==[run]== Funcotator CANONICAL con re-entrada inteligente (salida: ${OUT_DIR})"
for V in "${VCF_LIST[@]}"; do
  base_no_gz="${V%.vcf.gz}"
  bn="$(basename "${base_no_gz}")"
  out_vcf="${OUT_DIR}/${bn}.funcotator.vcf.gz"
  out_maf="${OUT_DIR}/${bn}.funcotator.maf"

  # Si ya existe el MAF, saltamos todo
  if [ -s "${out_maf}" ]; then
    echo "[skip] MAF ya existe: $(basename "${out_maf}")"
    continue
  fi

  # Si NO existe el VCF anotado, lo generamos
  if [ ! -s "${out_vcf}" ]; then
    echo "[VCF]  $(basename "${V}") -> $(basename "${out_vcf}")"
    gatk --java-options "-Xmx${JAVAMEM_GB}g" Funcotator \
      --variant "${V}" \
      --reference "${REF_FA}" \
      --ref-version "${REFVER}" \
      --data-sources-path "${DS_BASE}" \
      --output "${out_vcf}" \
      --output-file-format VCF \
      --transcript-selection-mode CANONICAL \
      --verbosity WARNING
    tabix -f "${out_vcf}" || true
  else
    echo "[skip] VCF anotado ya existe: $(basename "${out_vcf}")"
  fi

  # Generar MAF (si falta)
  if [ ! -s "${out_maf}" ]; then
    echo "[MAF]  $(basename "${V}") -> $(basename "${out_maf}")"
    gatk --java-options "-Xmx${JAVAMEM_GB}g" Funcotator \
      --variant "${V}" \
      --reference "${REF_FA}" \
      --ref-version "${REFVER}" \
      --data-sources-path "${DS_BASE}" \
      --output "${out_maf}" \
      --output-file-format MAF \
      --transcript-selection-mode CANONICAL \
      --annotation-default Center:OrganoidLab \
      --verbosity WARNING
  fi
done

echo "==[done]== Resultados en ${OUT_DIR}:"
ls -lh "${OUT_DIR}" | sed -n "1,200p"
'@

Invoke-DockerScript -Volume $VolumeName -Image 'broadinstitute/gatk:4.5.0.0' -ScriptBody $annotSh -Env @{ REFVER = $RefVersion; JAVAMEM_GB = $JavaGB }

# Copia al host si se pidió
if ($ExportDir) {
  $copySh = @'
set -e
mkdir -p /hostout
cp -av /work/funcotator_out/. /hostout/
echo "== Copiado a /hostout =="
ls -lh /hostout | sed -n "1,200p"
'@
  Invoke-DockerScript -Volume $VolumeName -Image 'alpine' -ScriptBody $copySh -ExtraMount "$ExportDir`:/hostout" -CPUs_ 1 -MemGB_ 1
  Write-Host "==> Copia completada en: $ExportDir"
}

Write-Host "Listo."
