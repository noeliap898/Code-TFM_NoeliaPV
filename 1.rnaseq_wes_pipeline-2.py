# ----------------------------------------------------------------------------------
# Paso 1: Clasificación automática de muestras
# ----------------------------------------------------------------------------------

import os
from collections import defaultdict

input_dir = "/Volumes/SecTFM"
muestras = defaultdict(lambda: defaultdict(dict))   


# Recorre el directorio de entrada y clasifica los FASTQ por paciente/clase/pareja
for archivo in os.listdir(input_dir):
    if archivo.startswith("._"):
        continue  # Saltar archivos ocultos innecesarios (archivos fantasma de macOS)

    # Indicamos que solo procese ficheros FASTQ/FASTQ.gz 
    if archivo.endswith(".FASTQ") or archivo.endswith(".FASTQ.gz"):
        partes = archivo.split("~")                       # Convención de nombres separada por "~"
        paciente = partes[0]                              # ID paciente/órgano antes del primer "~"
        tipo = partes[-1].split(".")[0]                   # Último bloque antes del primer punto (p.ej., "WES" o "RNASEQ")
        pareja = partes[-1].split(".")[1]  # R1 o R2     # Segundo token tras el último "~", separado por "." (esperado: R1/R2)

        # Determina la clase de datos por substring en el nombre del archivo
        if "germlineWES" in archivo:
            clase = "germlineWES"       # WES germinal
        elif "RNASEQ" in archivo:
            clase = "rnaseq"            # RNA-seq
        elif "WES" in archivo:
            clase = "wes"               # WES tumor/organoide
        else:
            continue                    # Si no contiene ninguna clase conocida, se ignora

        # Guarda la ruta absoluta del archivo en la estructura paciente->clase->(R1/R2)
        muestras[paciente][clase][pareja] = os.path.join(input_dir, archivo)
        
# Limpiar entradas no válidas que empiezan por "._"
muestras_limpias = {k: v for k, v in muestras.items() if not k.startswith("._")}
        
        
print(muestras_limpias)         #Así ya tenemos las muestras clasificadas segun sean de RNASeq, WES germinal o WES del organoide



# ----------------------------------------------------------------------------------
# Paso 2: Control de calidad de datos de entrada
# Lanza FastQC sobre todos los FASTQ y luego agrega resultados con MultiQC
# ----------------------------------------------------------------------------------

import subprocess

input_dir = "/Volumes/SecTFM"                               
qc_dir = "/Volumes/SecTFM/resultados/calidad_fastq"         # Carpeta de salida para QC
os.makedirs(qc_dir, exist_ok=True)

# Buscar todos los archivos FASTQ (exactamente con estas extensiones y mayúsculas)
fastq_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".FASTQ") or f.endswith(".FASTQ.gz")]

# Ejecutar FastQC en todos los FASTQ encontrados
for f in fastq_files:
    subprocess.run(["fastqc", f, "-o", qc_dir])

# Agrupar todos los informes con MultiQC sobre la carpeta de QC
subprocess.run(["multiqc", qc_dir, "-o", qc_dir])


#¿Qué se está generando?
##Por cada archivo .FASTQ.gz, se generan dos archivos:
###nombre_muestra_fastqc.html: informe gráfico con métricas como calidad de bases, contenido GC, adaptadores, etc.
###nombre_muestra_fastqc.zip: contiene los datos sin procesar y gráficos
##Y luego, MultiQC genera:
###multiqc_report.html: resumen agrupando todas las muestras en una sola vista
###Carpetas con datos agregados (multiqc_data/)

