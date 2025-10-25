

"""
Integración mutación + expresión (doble fila por gen TNBC con separador gris)
- Mutación (0/1) -> blanco / rojo
- Expresión (z-score log1p TPM) -> RdBu_r
- Separador gris
- Leyendas colocadas en vertical (una arriba de otra, sin solaparse)
"""

import argparse, os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

DEFAULT_TNBC = [
    "TP53","BRCA1","BRCA2","RB1","PTEN","PIK3CA","CDKN2A","ARID1A","ATM","ATR",
    "PALB2","FANCD2","RAD51","BARD1","CHEK1","CHEK2","NF1","KMT2C","KMT2D",
    "NOTCH1","NOTCH2","NOTCH3","EGFR","ERBB2","ERBB3","AKT1","AKT2","AKT3",
    "MTOR","TSC1","TSC2","MAP3K1","MAP2K4","KRAS","NRAS","MYC","NBN","PPM1D"
]

def pick_gene_col(df):
    for c in ["gene","Gene","Hugo_Symbol"]:
        if c in df.columns:
            return c
    return df.columns[0]

def simplify_sample_id(col: str) -> str:
    m = re.match(r"^(\d{3,})", col)
    if m:
        return m.group(1)
    return col.split("~")[0]

def zscore_by_row(M: pd.DataFrame) -> pd.DataFrame:
    X = M.values.astype(float)
    mu = np.nanmean(X, axis=1, keepdims=True)
    sd = np.nanstd(X, axis=1, ddof=0, keepdims=True)
    sd[sd==0] = 1.0
    Z = (X - mu) / sd
    return pd.DataFrame(Z, index=M.index, columns=M.columns)

def coerce_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out

def main(tpm_tsv, mut_tsv, out_png, tnbc_panel=None, zcap=2.0, sep_color="#4d4d4d"):
    # leer TPM
    tpm = pd.read_csv(tpm_tsv, sep="\t")
    gcol_tpm = pick_gene_col(tpm)
    tpm[gcol_tpm] = tpm[gcol_tpm].astype(str)
    tpm = tpm.drop_duplicates(subset=[gcol_tpm]).set_index(gcol_tpm)

    # leer MUT
    mut = pd.read_csv(mut_tsv, sep="\t")
    gcol_mut = pick_gene_col(mut)
    mut[gcol_mut] = mut[gcol_mut].astype(str)
    mut = mut.drop_duplicates(subset=[gcol_mut]).set_index(gcol_mut)

    # simplificar IDs
    tpm.columns = [simplify_sample_id(c) for c in tpm.columns]
    mut.columns = [simplify_sample_id(c) for c in mut.columns]

    # intersección
    common_samples = [s for s in tpm.columns if s in mut.columns]
    if len(common_samples) < 2:
        raise SystemExit(f"[ERROR] No hay suficientes muestras comunes: {common_samples}")
    tpm = tpm[common_samples]
    mut = mut[common_samples]

    # panel TNBC
    if tnbc_panel and os.path.exists(tnbc_panel):
        with open(tnbc_panel,"r",encoding="utf-8") as fh:
            genes_tnbc = [g.strip() for g in fh if g.strip()]
    else:
        genes_tnbc = list(DEFAULT_TNBC)

    genes_final = [g for g in genes_tnbc if g in tpm.index and g in mut.index]
    if not genes_final:
        raise SystemExit("[ERROR] Ninguno de los genes TNBC está en ambas matrices.")

    # preparar matrices numéricas
    tpm_num = coerce_numeric_df(tpm.loc[genes_final])
    expr_z = zscore_by_row(np.log1p(tpm_num))
    mut_num = coerce_numeric_df(mut.loc[genes_final]).clip(0,1)

    # construir bloques (mut, expr, separador) por gen
    blocks = []
    labels = []
    for g in genes_final:
        blocks.extend([mut_num.loc[g].values, expr_z.loc[g].values, np.zeros(len(common_samples)) * np.nan])
        labels.extend([f"{g} (mut)", f"{g} (expr)", ""])

    M = np.vstack(blocks)

    # figura
    plt.figure(figsize=(12, max(6, 0.55*len(genes_final))))
    ax = plt.gca()
    divider = make_axes_locatable(ax)

    # capas: separador gris, expresión, mutación
    sep_mask = np.full_like(M, np.nan, dtype=float)
    sep_mask[2::3, :] = 0.0
    cmap_sep = ListedColormap([sep_color])
    ax.imshow(sep_mask, aspect="auto", cmap=cmap_sep, vmin=0, vmax=1, interpolation="none")

    expr_mask = np.full_like(M, np.nan, dtype=float)
    expr_mask[1::3, :] = M[1::3, :]
    im_expr = ax.imshow(expr_mask, aspect="auto", cmap="RdBu_r", vmin=-zcap, vmax=zcap, interpolation="none")

    mut_mask = np.full_like(M, np.nan, dtype=float)
    mut_mask[0::3, :] = M[0::3, :]
    cmap_mut = ListedColormap(["white","#b22222"])
    im_mut = ax.imshow(mut_mask, aspect="auto", cmap=cmap_mut, vmin=0, vmax=1, interpolation="none")

    # ejes y etiquetas
    ax.set_xticks(np.arange(len(common_samples)))
    ax.set_xticklabels(common_samples, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("ID Muestra")
    ax.set_ylabel("Genes (mutación / expresión)")
    ax.set_title("Integración mutación y expresión en genes TNBC")

    # leyendas apiladas
    cax_expr = divider.append_axes("right", size="3%", pad=0.3)
    cb1 = plt.colorbar(im_expr, cax=cax_expr)
    cb1.set_label("z-score (log1p TPM)")

    cax_mut = divider.append_axes("right", size="3%", pad=1.0)  # más arriba
    cb2 = plt.colorbar(im_mut, cax=cax_mut)
    cb2.set_label("Mutación (0=no, 1=sí)")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print(f"[OK] Guardado -> {out_png}")

if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpm", required=True)
    ap.add_argument("--mut", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tnbc-panel", default=None)
    ap.add_argument("--zcap", type=float, default=2.0)
    args = ap.parse_args()
    main(args.tpm, args.mut, args.out, args.tnbc_panel, args.zcap)
