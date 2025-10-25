# ----------------------------------------------------------------------------------
# Paso 4: MultiQC post limpieza
# ----------------------------------------------------------------------------------
###No renombra nada
###En fastqc_pre/ puedes tener .html y .zip; el script solo copia los _fastqc.zip
###Copia los PRE a qc_all/PRE/ y los POST a qc_all/POST/
###Lanza MultiQC sobre qc_all/ con --dirs --dirs-depth 1 --fullnames para que el informe muestre claramente PRE/ y POST/ y los nombres completos de archivo (p. ej., …FASTQ_fastqc.zip vs …FASTQ.gz_val_1_fastqc.zip)



import os
import shutil
import subprocess

BASE = "/Volumes/SecTFM"
PRE_DIR  = os.path.join(BASE, "fastqc_pre")      # aquí tienes .html y .zip (PRE)
POST_DIR = os.path.join(BASE, "cleaned_reads")   # aquí están los .zip POST (TrimGalore)
QC_ALL   = os.path.join(BASE, "qc_all")
QC_PRE   = os.path.join(QC_ALL, "PRE")
QC_POST  = os.path.join(QC_ALL, "POST")

def is_fastqc_zip(fn: str) -> bool:
    return fn.endswith("_fastqc.zip")

def copy_fastqc_zips(src_dir: str, dst_dir: str) -> int:
    os.makedirs(dst_dir, exist_ok=True)
    n = 0
    for fn in os.listdir(src_dir):
        if is_fastqc_zip(fn):
            shutil.copy2(os.path.join(src_dir, fn), os.path.join(dst_dir, fn))
            n += 1
    return n

##Copiar SOLO los .zip de FastQC 
if not os.path.isdir(PRE_DIR):
    raise SystemExit(f"Directorio PRE no encontrado: {PRE_DIR}")
if not os.path.isdir(POST_DIR):
    raise SystemExit(f"Directorio POST no encontrado: {POST_DIR}")

os.makedirs(QC_ALL, exist_ok=True)

n_pre  = copy_fastqc_zips(PRE_DIR,  QC_PRE)   # copia solo *_fastqc.zip (ignora .html)
n_post = copy_fastqc_zips(POST_DIR, QC_POST)

print(f"Copiados PRE  (.zip) -> {QC_PRE}:  {n_pre}")
print(f"Copiados POST (.zip) -> {QC_POST}: {n_post}")

##Revisar longitudes
count_qc_pre  = len([f for f in os.listdir(QC_PRE)  if is_fastqc_zip(f)])
count_qc_post = len([f for f in os.listdir(QC_POST) if is_fastqc_zip(f)])
print(f"Verificación: PRE={count_qc_pre} | POST={count_qc_post}")

##Ejecutar MultiQC mostrando nombres completos y carpeta padre (PRE/POST)
# Requiere: pip install multiqc
multiqc_cmd = [
    "multiqc",
    QC_ALL,                         # analiza qc_all con PRE/ y POST/
    "-o", QC_ALL,                   # deja el informe dentro de qc_all
    "--filename", "multiqc_pre_post.html",
    "--title", "QC PRE vs POST (FastQC)",
    "--dirs",                       # añade el nombre del directorio a la muestra (PRE/POST)
    "--dirs-depth", "1",
    "--fullnames"                   # usa el nombre completo del archivo (no lo “limpia”)
]

print("Lanzando MultiQC...")
res = subprocess.run(multiqc_cmd)
if res.returncode == 0:
    print(f"Informe generado: {os.path.join(QC_ALL, 'multiqc_pre_post.html')}")
else:
    print("Error al ejecutar MultiQC. Asegúrate de tenerlo instalado: `pip install multiqc`")
