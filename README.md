#Validación de organoides derivados de pacientes como modelos representativos en cáncer de mama triple negativo mediante análisis multi-ómico

**Autora:** Noelia Pérez Ventura  
**Director:** Francisco Javier Ruiz Ojeda  
**Universidad Internacional de La Rioja (UNIR)**  
**Máster Universitario en Bioinformática – Septiembre 2025**

---

## Resumen

El cáncer de mama triple negativo (CMTN) constituye uno de los subtipos más agresivos y heterogéneos de esta enfermedad, caracterizado por la ausencia de receptores hormonales y HER2, lo que limita las opciones terapéuticas y empeora el pronóstico.  
El presente trabajo tuvo como objetivo **validar organoides derivados de pacientes (PDOs)** como modelos representativos del CMTN mediante un enfoque **multi-ómico**, combinando datos genómicos (WES) y transcriptómicos (RNA-seq).

Se analizaron **seis pares organoide–paciente** del *Patient-Derived Models Repository (PDMR)* del *National Cancer Institute (NCI)*, con WES germinal, WES tumoral del organoide y RNA-seq del mismo.  
El pipeline permitió comparar la carga mutacional, el perfil de expresión génica y la concordancia molecular entre organoides y tumores originales, demostrando una alta estabilidad genómica y transcriptómica de los PDOs.

Los resultados confirman que los organoides constituyen modelos preclínicos robustos y versátiles, útiles para estudios funcionales, cribado de fármacos y estrategias de medicina personalizada en cáncer de mama triple negativo.

---

## Estructura del repositorio

El repositorio contiene los scripts empleados en el flujo completo de análisis multi-ómico, agrupados por etapas:

### 1. Preprocesamiento y control de calidad (WES y RNA-seq)

| Script | Descripción |
|--------|--------------|
| `1.rnaseq_wes_pipeline-2.py` | Define el flujo general de análisis WES + RNA-seq, orquestando la ejecución de los módulos subsiguientes. |
| `2.1. Limpieza de secuencias post FASTQC.py` | Limpieza de lecturas crudas de WES y RNA-seq con TrimGalore/Cutadapt tras control de calidad inicial con FastQC. |
| `2.2. MultiQC post limpieza.py` | Generación de informes agregados de calidad tras el recorte, mediante MultiQC. |

---

### 2. Procesamiento WES (Whole-Exome Sequencing)

| Script | Descripción |
|--------|--------------|
| `3.prepare_wesdata.ps1` | Organización de directorios, verificación de FASTQ y preparación de inputs para alineamiento. |
| `4.wes_align_all_bwa_fix4.ps1` | Alineamiento de lecturas WES al genoma humano GRCh38 usando BWA-MEM2. |
| `5.wes_markduplicates.ps1` | Marcado y eliminación de duplicados con GATK MarkDuplicatesSpark. |
| `6.wes_somatic_mutect2_dedup.ps1` | Llamada de variantes somáticas con GATK Mutect2, usando muestras germinales y tumorales pareadas. |
| `7.wes_filter_mutect.ps1` | Filtrado de variantes somáticas con GATK FilterMutectCalls. |
| `8.wes_mutect2_qc.ps1` | Evaluación de métricas de calidad de variantes y generación de informes QC. |
| `9.wes_funcotator_download_and_run.ps1` | Descarga de recursos y anotación funcional de variantes con GATK Funcotator. |
| `10.run_funcotator_filtered_canonical.ps1` | Filtrado final de variantes anotadas y exportación en formato MAF para análisis posterior. |
| `11.build_mutation_matrix_softfilters_v3.py` | Construcción de matrices mutacionales por muestra y gen, integrando filtros suaves. |
| `12.wes_filtering_qc_v2.py` | Control de calidad adicional y generación de resúmenes de carga mutacional (TMB). |

---

### 3. Procesamiento RNA-seq

| Script | Descripción |
|--------|--------------|
| `13.run_rnaseq_salmon_pairs.ps1` | Cuantificación de la expresión génica en modo paired-end con **Salmon**, obteniendo matrices TPM. |
| `14.harmonize_ids_and_save.py` | Homogeneización de identificadores génicos (Ensembl ↔ GENCODE) y consolidación de matrices de expresión. |

---

### 4. Integración multi-ómica y visualización

| Script | Descripción |
|--------|--------------|
| `Fig 2_fig_heatmap_mutations_variable.py` | Generación del *heatmap* de mutaciones somáticas más variables (Figura 2 del TFM). |
| `fig3_expr_heatmap_TNBCfirst.py` | *Heatmap* de genes de expresión diferencial entre organoides (Figura 3). |
| `fig4_top_genes_per_sample_v2.py` | Visualización de los genes más expresados por muestra (Figura 4). |
| `fig5_pca_expression_v2.py` | Análisis de componentes principales (PCA) de expresión génica (Figura 5). |
| `fig6_sample_correlation_matrix_clean.py` | Cálculo y representación de la matriz de correlación entre organoides (Figura 6). |
| `fig7_integracion_TNBC_doble_spaced_v2.py` | Integración de mutaciones y expresión para genes característicos del CMTN (Figura 7). |

---

## Requisitos del entorno

El proyecto se desarrolló y ejecutó en entornos **Docker** y **Conda** para garantizar la reproducibilidad.

### Dependencias principales

- **Python** ≥ 3.10  
- **PowerShell** ≥ 7.2  
- **Docker** ≥ 24  
- **Conda (Miniconda/Anaconda)** ≥ 23  
- **GATK** 4.4  
- **BWA-MEM2** 2.2.1  
- **FastQC** 0.12.1  
- **TrimGalore / Cutadapt** 5.1  
- **MultiQC** 1.15  
- **Salmon** 1.10  
- **Funcotator** (bundles GENCODE v43, dbSNP 151, ClinVar)  
- **Pandas**, **NumPy**, **Matplotlib**, **Seaborn**, **scikit-learn**

---

##Instalación y ejecución

### 1️. Clonar el repositorio

```bash
git clone https://github.com/tu_usuario/tu_repositorio.git
cd tu_repositorio
```

### 2️. Crear el entorno Conda

```bash
conda create -n tn_bc_multiomics python=3.10
conda activate tn_bc_multiomics
pip install pandas numpy matplotlib seaborn scikit-learn multiqc
```

### 3️. Configurar Docker (para GATK, BWA, Funcotator)

```bash
docker pull broadinstitute/gatk:4.4.0.0
docker pull combinelab/salmon:1.10.0
```

### 4️. Ejecutar el flujo completo

1. **Preprocesamiento WES y RNA-seq**
   ```bash
   python 2.1.\ Limpieza\ de\ secuencias\ post\ FASTQC.py
   python 2.2.\ MultiQC\ post\ limpieza.py
   ```
2. **Alineamiento y llamada de variantes**
   ```powershell
   pwsh ./4.wes_align_all_bwa_fix4.ps1
   pwsh ./6.wes_somatic_mutect2_dedup.ps1
   pwsh ./9.wes_funcotator_download_and_run.ps1
   ```
3. **RNA-seq**
   ```powershell
   pwsh ./13.run_rnaseq_salmon_pairs.ps1
   python 14.harmonize_ids_and_save.py
   ```
4. **Integración multi-ómica y visualización**
   ```bash
   python Fig\ 2_fig_heatmap_mutations_variable.py
   python fig7_integracion_TNBC_doble_spaced_v2.py
   ```

Los resultados (VCF, MAF, TPM, figuras y métricas QC) se guardan automáticamente en los subdirectorios `/results` y `/figures`.

---

## Resultados esperados

- **Carga mutacional somática (TMB)** por organoide y comparación con sus controles germinales.  
- **Perfiles transcriptómicos** expresados como TPM, con *heatmaps*, PCA y correlaciones.  
- **Integración mutación–expresión** mostrando la conservación funcional de vías clave en el CMTN (PI3K/AKT, TP53, BRCA1).  
- **Visualizaciones reproducibles** correspondientes a las Figuras 2–7 del TFM.

---

## Autoría y agradecimientos

Este repositorio acompaña al Trabajo Fin de Máster:

> **“Validación de organoides derivados de pacientes como modelos representativos en cáncer de mama triple negativo mediante análisis multi-ómico”**  
> Presentado por **Noelia Pérez Ventura**  
> Máster Universitario en Bioinformática – Universidad Internacional de La Rioja (UNIR)  
> Septiembre de 2025  

Agradecimientos al **Prof. Francisco Javier Ruiz Ojeda** por la dirección y supervisión del proyecto, y al **National Cancer Institute (PDMR)** por la disponibilidad de los datos de secuenciación.

---
