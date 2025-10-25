

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def sample_id_only(colname: str) -> str:
    """
    Extrae el ID paciente (primer bloque antes del primer '~').
    Si ya está limpio, lo deja tal cual.
    """
    if "~" in colname:
        return colname.split("~", 1)[0]
    return colname

def load_tpm_matrix(tpm_tsv: str) -> pd.DataFrame:
    M = pd.read_csv(tpm_tsv, sep="\t", index_col=0)
    # A numérico (si hay NA/strings raros -> 0)
    M = M.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # Renombra columnas a solo ID
    new_cols = [sample_id_only(c) for c in M.columns]
    M.columns = new_cols
    # Si hay columnas duplicadas tras “limpiar” nombres, haz la media
    if len(set(new_cols)) != len(new_cols):
        M = M.groupby(level=0, axis=1).mean()
    return M

def top_genes_per_sample(M: pd.DataFrame, topk: int, log1p: bool) -> dict:
    """
    Devuelve un dict {sample_id: DataFrame con columnas [gene, value]}
    con los top-k genes por muestra.
    """
    results = {}
    for s in M.columns:
        vals = M[s].copy()
        if log1p:
            vals = np.log1p(vals)
        top = vals.sort_values(ascending=False).head(topk)
        df = top.reset_index()
        df.columns = ["gene", "value"]
        results[s] = df
    return results

def plot_panels(tops: dict, out_png: str, topk: int, log1p: bool):
    samples = sorted(tops.keys())
    n = len(samples)
    # cuadrícula 2x3 (sirve para tus 6 organoides); si cambia n, ajusta automáticamente
    ncols = 3
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(15, 9),
        constrained_layout=True
    )
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])

    # coloreado simple coherente
    for i, s in enumerate(samples):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        df = tops[s]
        # barras horizontales
        ax.barh(df["gene"], df["value"])
        ax.invert_yaxis()  # mayor arriba
        ax.set_title(s, fontsize=12, fontweight="bold")
        ax.set_xlabel("log1p(TPM)" if log1p else "TPM")
        ax.set_ylabel("Genes")
        # tipografías legibles
        ax.tick_params(axis='y', labelsize=8)
        ax.tick_params(axis='x', labelsize=9)

    # oculta ejes vacíos si n no llena la rejilla
    for j in range(n, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    suptitle = f"Genes más expresados por organoide"
    if log1p:
        suptitle += " (escala log1p)"
    fig.suptitle(suptitle, fontsize=14, y=1.02)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[OK] Figura 6 generada -> {out_png}")

def main():
    ap = argparse.ArgumentParser(
        description="Figura 6: Top genes más expresados por organoide (barras horizontales)."
    )
    ap.add_argument("--tpm", required=True, help="Matriz TPM (genes en filas, muestras en columnas)")
    ap.add_argument("--out", required=True, help="PNG de salida")
    ap.add_argument("--topk", type=int, default=12, help="Nº de genes por muestra (default: 12)")
    ap.add_argument("--log1p", action="store_true", help="Usar log1p(TPM) para estabilizar rangos")
    args = ap.parse_args()

    M = load_tpm_matrix(args.tpm)
    tops = top_genes_per_sample(M, topk=args.topk, log1p=args.log1p)
    plot_panels(tops, args.out, topk=args.topk, log1p=args.log1p)

if __name__ == "__main__":
    main()
