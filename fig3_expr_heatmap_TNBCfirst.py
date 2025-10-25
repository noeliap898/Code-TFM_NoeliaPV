

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# TNBC por defecto (mismo panel que en mutaciones)
DEFAULT_TNBC = [
    "TP53","BRCA1","BRCA2","RB1","CDKN2A","PTEN","PIK3CA","ARID1A",
    "ATR","PALB2","ATM","CHEK1","CHEK2","RAD51","RAD51C","RAD51D","BARD1",
    "NBN","FANCA","FANCD2","FANCI","FANCM","RAD50","MRE11","BLM",
    "EGFR","ERBB2","ERBB3","AKT1","AKT2","AKT3","mTOR","TSC1","TSC2",
    "MAP2K4","MAP3K1","KRAS","NRAS","MYC","KMT2C","KMT2D",
    "NOTCH1","NOTCH2","NOTCH3"
]

def simplify_sample_id(x: str) -> str:
    """
    Extrae el ID numérico inicial (antes del primer '~') si existe.
    Si no, devuelve el string original.
    """
    return x.split("~", 1)[0]

def load_tnbc_order(order_file: str|None):
    if order_file is None:
        return DEFAULT_TNBC[:]
    # se espera un archivo con un gen por línea
    genes = []
    with open(order_file, "r", encoding="utf-8") as fh:
        for line in fh:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.append(g)
    # si el archivo existe pero está vacío, caemos al panel por defecto
    return genes if genes else DEFAULT_TNBC[:]

def zscore_by_row(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=1)
    sd = df.std(axis=1).replace(0, np.nan)
    return (df.sub(mu, axis=0)).div(sd, axis=0)

def main(tpm_tsv, out_png, topn=150, tnbc_order_file=None):
    # carga TPM
    df = pd.read_csv(tpm_tsv, sep="\t", index_col=0)
    # a numérico silencioso
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # renombra columnas a solo ID (388244, etc.)
    df.columns = [simplify_sample_id(c) for c in df.columns]

    # panel TNBC (en orden)
    tnbc_list = load_tnbc_order(tnbc_order_file)
    # quedarnos solo con genes presentes
    tnbc_present = [g for g in tnbc_list if g in df.index]

    # si faltan muchos genes, avisamos en consola
    missing = sorted(set(tnbc_list) - set(tnbc_present))
    if missing:
        print(f"[WARN] {len(missing)} genes TNBC no están en la matriz TPM (se omiten). Ejemplos: {missing[:8]}")

    # submatriz TNBC
    M_tnbc = df.loc[tnbc_present] if len(tnbc_present) else pd.DataFrame(index=[], columns=df.columns)

    # genes más variables (excluyendo TNBC)
    resto = df.drop(index=tnbc_present, errors="ignore")
    var = resto.var(axis=1)
    # quita los que tienen var = 0 (constantes)
    var = var[var > 0]
    sel = var.sort_values(ascending=False).head(topn).index
    M_var = resto.loc[sel]

    # z-score por gen (fila) para cada panel 
    Mz_tnbc = zscore_by_row(M_tnbc) if len(M_tnbc) else M_tnbc
    Mz_var  = zscore_by_row(M_var)  if len(M_var)  else M_var

    # figura: dos paneles apilados, con ejes y etiquetas 
    nrows_top = max(2, len(Mz_tnbc))
    nrows_bot = max(2, len(Mz_var))
    # altura proporcional a nº de genes para que se lean todas las etiquetas
    h_top = 0.28 * nrows_top
    h_bot = 0.28 * nrows_bot
    fig_h = min(45, h_top + h_bot + 2.5)

    fig, axes = plt.subplots(
        nrows=2, ncols=1, figsize=(10, fig_h),
        gridspec_kw={"height_ratios": [h_top, h_bot]}
    )

    # panel superior
    ax1 = axes[0]
    if len(Mz_tnbc):
        sns.heatmap(
            Mz_tnbc, cmap="vlag", center=0, cbar_kws={"label": "z-score"},
            ax=ax1
        )
    ax1.set_title("Panel de expresión de genes TNBC")
    ax1.set_xlabel("")  # sin etiqueta en el superior
    ax1.set_ylabel("Genes")

    # panel inferior
    ax2 = axes[1]
    if len(Mz_var):
        sns.heatmap(
            Mz_var, cmap="vlag", center=0, cbar_kws={"label": "z-score"},
            ax=ax2
        )
    ax2.set_title("Panel de expresión de genes con mayor variabilidad")
    ax2.set_xlabel("ID Muestra")
    ax2.set_ylabel("Genes")

    # rotación de ticks para que quepan
    for ax in axes:
        ax.tick_params(axis='x', labelrotation=45, labelright=False)
        ax.tick_params(axis='y', labelrotation=0)

    # línea divisoria visible entre paneles (borde grueso del superior)
    for spine in ["bottom"]:
        ax1.spines[spine].set_visible(True)
        ax1.spines[spine].set_linewidth(3)

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300)
    print(f"[OK] Heatmap de expresión guardado -> {out_png}")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Heatmap de expresión: TNBC (arriba) + genes más variables (abajo)")
    ap.add_argument("--tpm", required=True, help="Matriz de TPM por gen (filas) × muestra (columnas) en TSV")
    ap.add_argument("--out_png", required=True, help="Ruta de salida PNG")
    ap.add_argument("--topn", type=int, default=150, help="Nº de genes variables (excluyendo TNBC) para el panel inferior")
    ap.add_argument("--tnbc_order", default=None, help="(Opcional) fichero de texto con el orden TNBC usado en mutaciones (1 gen/linea)")
    args = ap.parse_args()
    main(args.tpm, args.out_png, args.topn, args.tnbc_order)
