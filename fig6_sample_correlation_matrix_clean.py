

import argparse, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def clean_id(x: str) -> str:
    """
    Convierte nombres tipo:
      388244~064-R~V1-organoid~v2.0.2.51.0
    en:
      388244
    """
    if pd.isna(x):
        return x
    # Toma lo que haya antes del primer '~' si existe, si no, deja el número/ID
    m = re.match(r"^([0-9]+)", str(x))
    if m:
        return m.group(1)
    # fallback: corta por '~' y toma el primer trozo
    return str(x).split("~")[0]

def load_tpm(path: str) -> pd.DataFrame:
    # lee matriz gen x muestra con primera columna como gene_name si existe
    df = pd.read_csv(path, sep=None, engine="python")
    # detectar columna de genes
    gene_col = None
    for c in ["gene", "gene_id", "gene_name", "Gene", "GeneSymbol", "symbol"]:
        if c in df.columns:
            gene_col = c
            break
    if gene_col is None:
        # si no se encontró, asume la primera columna
        gene_col = df.columns[0]
    df = df.set_index(gene_col)
    # eliminar posibles filas vacías o duplicadas
    df = df[~df.index.isna()].copy()
    df = df[~df.index.duplicated(keep="first")].copy()
    return df

def make_correlation(df_tpm: pd.DataFrame,
                     out_png: str,
                     out_csv: str,
                     annotate: bool = False,
                     cluster: bool = False):
    # limpiar nombres de muestras
    original_cols = df_tpm.columns.tolist()
    cleaned_cols = [clean_id(c) for c in original_cols]
    # si hay colisiones de ID (dos columnas acaban con el mismo ID), añade sufijo
    counts = {}
    final_cols = []
    for c in cleaned_cols:
        if c not in counts:
            counts[c] = 1
            final_cols.append(c)
        else:
            counts[c] += 1
            final_cols.append(f"{c}_{counts[c]}")
    df_tpm.columns = final_cols

    # quitar genes con todo cero (evita correlaciones indefinidas)
    M = df_tpm.copy()
    M = M.loc[(M.sum(axis=1) > 0), :]

    # log1p para estabilizar (opcional, útil para TPM)
    M_log = np.log1p(M)

    # correlación (Pearson)
    corr = M_log.corr(method="pearson")

    # guardar CSV con IDs limpios
    corr.to_csv(out_csv, sep="\t", float_format="%.3f")

    # figura
    plt.figure(figsize=(8, 6))
    cmap = sns.diverging_palette(240, 10, as_cmap=True)  # azul-blanco-rojo
    # opcional clustering
    if cluster:
        g = sns.clustermap(
            corr,
            cmap=cmap,
            vmin=0, vmax=1,
            linewidths=0.5, linecolor="white",
            annot=annotate, fmt=".2f",
            cbar_kws={"label": "Correlación de Pearson"},
            figsize=(8, 8)
        )
        # títulos y ejes
        g.ax_heatmap.set_title("Correlación entre organoides (TPM, log1p)", pad=16)
        g.ax_heatmap.set_xlabel("ID Muestra")
        g.ax_heatmap.set_ylabel("ID Muestra")
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        ax = sns.heatmap(
            corr,
            cmap=cmap,
            vmin=0, vmax=1,
            square=True,
            linewidths=0.5, linecolor="white",
            annot=annotate, fmt=".2f",
            cbar_kws={"label": "Correlación de Pearson"},
        )
        ax.set_title("Correlación entre organoides (TPM, log1p)", pad=12)
        ax.set_xlabel("ID Muestra")
        ax.set_ylabel("ID Muestra")
        plt.xticks(rotation=0)
        plt.yticks(rotation=0)
        plt.tight_layout()
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()

    print(f"[OK] Guardado heatmap: {out_png}")
    print(f"[OK] Guardada matriz TSV: {out_csv}")

def main():
    ap = argparse.ArgumentParser(description="Matriz de correlaciones entre organoides (TPM).")
    ap.add_argument("--tpm", required=True, help="Ruta a tpm_gene_matrix.tsv (genes x muestras).")
    ap.add_argument("--out_png", required=True, help="Figura de salida (PNG).")
    ap.add_argument("--out_tsv", required=True, help="TSV de salida con la matriz de correlación.")
    ap.add_argument("--annot", action="store_true", help="Anotar los coeficientes en cada celda.")
    ap.add_argument("--cluster", action="store_true", help="Aplicar clustering jerárquico a las muestras.")
    args = ap.parse_args()

    tpm = load_tpm(args.tpm)
    make_correlation(tpm, args.out_png, args.out_tsv, annotate=args.annot, cluster=args.cluster)

if __name__ == "__main__":
    main()
