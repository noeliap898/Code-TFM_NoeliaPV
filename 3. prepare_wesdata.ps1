
param(
  [string]$VolumeName    = "wesdata",
  [string]$FuncotatorDir = "D:/refs/funcotator",   # Data sources (extraídos o .tar.gz)
  [string]$GenomeDir     = "D:/Datos",             # FASTA + .dict + .fai
  [string]$VcfDir        = "D:/results/wes"        # VCFs filtrados
)

Write-Host "Preparando volumen: $VolumeName"
Write-Host "DataSources (Windows): $FuncotatorDir"
Write-Host "Referencia  (Windows): $GenomeDir"
Write-Host "VCFs        (Windows): $VcfDir"

#Crear volumen si no existe
docker volume create $VolumeName | Out-Null

function Invoke-DockerScript {
  param(
    [string]$Volume,
    [string]$HostMap,     # puede ir vacío
    [string]$ScriptBody
  )
  # Normaliza a LF y guarda en tmp .sh
  $bodyLF = $ScriptBody -replace "`r`n","`n"
  $tmpSh  = [System.IO.Path]::GetTempFileName() + ".sh"
  [System.IO.File]::WriteAllText($tmpSh, $bodyLF, [System.Text.Encoding]::ASCII)

  $args = @('run','--rm','-v',"$Volume`:/work",'-v',"$tmpSh`:/root/run.sh:ro")
  if ($HostMap) { $args += @('-v',"$HostMap`:/src") }
  $args += @('alpine','sh','/root/run.sh')

  $proc = Start-Process -FilePath 'docker' -ArgumentList $args -NoNewWindow -PassThru -Wait
  $code = $proc.ExitCode
  Remove-Item $tmpSh -ErrorAction SilentlyContinue
  if ($code -ne 0) { throw "docker run falló (ExitCode=$code)" }
}

#Copiar DATA SOURCES a /work/refs/funcotator
$shData = @'
mkdir -p /work/refs/funcotator
# Si hay un tar.gz, extraer; si no, copiar con tar streaming
if [ -f /src/funcotator_dataSources.*.tar.gz ]; then
  tar -xzf /src/funcotator_dataSources.*.tar.gz -C /work/refs/funcotator
else
  tar -C /src -cf - . 2>/dev/null | tar -C /work/refs/funcotator -xf - || true
fi
# Symlink canónico si existe una carpeta versionada hg38
latest=$(ls -d /work/refs/funcotator/funcotator_dataSources.v1.*.hg38.* 2>/dev/null | sort | tail -n1 || true)
if [ -n "$latest" ]; then ln -sfn "$latest" /work/refs/funcotator/funcotator_dataSources; fi
echo "== Datasources en /work/refs/funcotator =="
ls -l /work/refs/funcotator || true
[ -d /work/refs/funcotator/funcotator_dataSources/hg38 ] || [ -d /work/refs/funcotator/funcotator_dataSources/hg19 ] || echo "[WARN] No veo funcotator_dataSources/{hg38,hg19}"
'@
Invoke-DockerScript -Volume $VolumeName -HostMap "D:/refs/funcotator" -ScriptBody $shData

#Copiar GENOMA a /work/ref
$shRef = @'
mkdir -p /work/ref
tar -C /src -cf - . 2>/dev/null | tar -C /work/ref -xf -
echo "== Referencia en /work/ref =="
ls -l /work/ref || true
'@
Invoke-DockerScript -Volume $VolumeName -HostMap "D:/Datos" -ScriptBody $shRef

#Copiar VCFs (recursivo) a /work/m2
$shVcf = @'
mkdir -p /work/m2
tar -C /src -cf - . 2>/dev/null | tar -C /work/m2 -xf - || true
echo "== Muestra de VCFs .filtered en /work/m2 =="
find /work/m2 -type f -name "*.filtered.vcf.gz" | head -n 20 || true
count=$(find /work/m2 -type f -name "*.filtered.vcf.gz" | wc -l | tr -d " ")
echo "[INFO] Total *.filtered.vcf.gz encontrados: $count"
'@
Invoke-DockerScript -Volume $VolumeName -HostMap "D:/results/wes" -ScriptBody $shVcf

Write-Host "==> Volumen ${VolumeName}: copia finalizada."
Write-Host ">>> Verificando contenido dentro del volumen..."

#Verificación compacta
$shVerify = @'
echo "--- /work/refs/funcotator ---"
ls -l /work/refs/funcotator || true
echo "--- /work/refs/funcotator/funcotator_dataSources ---"
ls -l /work/refs/funcotator/funcotator_dataSources 2>/dev/null || true
echo "--- /work/ref ---"
ls -l /work/ref || true
echo "--- /work/m2 (primeros .filtered) ---"
find /work/m2 -maxdepth 3 -type f -name "*.filtered.vcf.gz" | head -n 20 || true
echo "--- Conteo VCFs .filtered ---"
find /work/m2 -type f -name "*.filtered.vcf.gz" | wc -l | tr -d " "
'@
Invoke-DockerScript -Volume $VolumeName -HostMap "" -ScriptBody $shVerify

