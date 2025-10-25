
import os, json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
from sklearn.decomposition import PCA

def make_sample_colors(samples, outdir):
    pal = sns.color_palette("tab10", n_colors=len(samples))
    mapping = {s: mcolors.to_hex(pal[i]) for i, s in enumerate(samples)}
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "sample_colors.json"), "w") as f:
        json.dump(mapping, f, indent=2)
    return mapping

def load_sample_colors(samples, outdir):
    pj = os.path.join(outdir, "sample_colors.json")
    if os.path.exists(pj):
        try:
            mp = json.load(open(pj))
            if all(s in mp for s in samples):
                return mp
        except Exception:
            pass
    return make_sample_colors(samples, outdir)

def simplify_ids(samples):
    """Extraer solo el ID numérico de cada muestra (ej. 388244)."""
    return [s.split("~")[0] for s in samples]

def main(tpm_tsv, out_png):
    df = pd.read_csv(tpm_tsv, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # log-transformación para estabilizar varianza
    X = np.log1p(df).T  # muestras x genes
    pca = PCA(n_components=2, random_state=1)
    pcs = pca.fit_transform(X)
    pc1, pc2 = pca.explained_variance_ratio_[0]*100, pca.explained_variance_ratio_[1]*100

    samples_full = X.index.tolist()
    samples = simplify_ids(samples_full)  # solo IDs

    colors = load_sample_colors(samples, os.path.dirname(out_png))
    c = [colors[s] for s in samples]

    plt.figure(figsize=(6,5))
    scatter = plt.scatter(pcs[:,0], pcs[:,1], s=100, c=c, edgecolor="k", linewidth=0.7)

    # leyenda con IDs
    for s in set(samples):
        plt.scatter([], [], c=colors[s], label=s, s=80, edgecolor="k")

    plt.xlabel(f"PC1 ({pc1:.1f}%)")
    plt.ylabel(f"PC2 ({pc2:.1f}%)")
    plt.title("PCA de perfiles de expresión (log1p TPM)")
    plt.legend(title="ID Muestra", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print("[OK] Figura PCA generada:", out_png)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpm", required=True, help="Matriz de TPM (genes x muestras)")
    ap.add_argument("--out", required=True, help="Ruta de salida PNG")
    args = ap.parse_args()
    main(args.tpm, args.out)
