
# ----------------------------------------------------------------------------------
# Paso 3: Ejecuta recorte de calidad con TrimGalore sobre lecturas paired-end y QC post-trimming
# ----------------------------------------------------------------------------------

##Parámetros
#Recorte por calidad: Phred < 20
#Descarte de lecturas cortas (< 35 pb)
#Ejecución de FastQC tras el trimming
#Compresión de resultados (.gz)
#Uso de 4 hilos



import os
import subprocess
from collections import defaultdict

input_dir = "/Volumes/SecTFM"  
output_dir = os.path.join(input_dir, "cleaned_reads")
os.makedirs(output_dir, exist_ok=True)  # crea la carpeta de salida si no existe

## Ruta al ejecutable de TrimGalore 
trim_galore_path = "/Users/noeliaperezventura/TrimGalore-0.6.10/trim_galore"

## Añadir igzip/pigz al PATH (Homebrew suele instalar en /opt/homebrew/bin)
env = os.environ.copy()
env["PATH"] = "/opt/homebrew/bin:" + env["PATH"]

## Detectar muestras emparejadas R1/R2
muestras = defaultdict(dict)

for archivo in os.listdir(input_dir):
    if not archivo.endswith((".FASTQ", ".FASTQ.gz")):
        continue  # salta todo lo que no sea FASTQ
    if "R1" in archivo:
        clave = archivo.replace("R1", "")  # clave común para emparejar R1/R2
        muestras[clave]["R1"] = os.path.join(input_dir, archivo)
    elif "R2" in archivo:
        clave = archivo.replace("R2", "")
        muestras[clave]["R2"] = os.path.join(input_dir, archivo)

## Ejecutar TrimGalore sobre cada par válido R1/R2
for clave, pares in muestras.items():
    if "R1" in pares and "R2" in pares:
        r1 = pares["R1"]
        r2 = pares["R2"]
        print(f"Procesando muestra: {clave}")

        cmd = [
            trim_galore_path,
            "--paired",       # indica que son lecturas paired-end
            "--quality", "20",# recorta bases con calidad <20
            "--length", "35", # descarta lecturas <35 pb
            "--fastqc",       # ejecuta FastQC tras el trimming
            "--cores", "4",   # usa 4 hilos
            "--gzip",         # comprime las salidas
            "-o", output_dir, # carpeta de salida
            r1, r2
        ]

        ## Lanza TrimGalore
        resultado = subprocess.run(cmd, env=env)
        if resultado.returncode != 0:
            print(f"Error al procesar {clave}")
    else:
        ## Si falta R1 o R2, se notifica y no se ejecuta
        print(f"Par incompleto para {clave}")


#¿Qué se está generando?
##Archivos FASTQ.gz recortados en carpeta "cleaned_reads/".
##Informes FastQC post-trimming (*.html, *.zip).