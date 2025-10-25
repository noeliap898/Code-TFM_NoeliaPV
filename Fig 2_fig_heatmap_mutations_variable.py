
import argparse, pandas as pd, seaborn as sns, matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--matrix", required=True, help="mutated_genes_matrix_full.tsv o .xls")
parser.add_argument("--out_png", required=True)
parser.add_argument("--topk", type=int, default=150, help="nº de genes más variables")
parser.add_argument("--tnbc_first", action="store_true", help="poner genes panel TNBC arriba")
args = parser.parse_args()

# lee matriz
if args.matrix.lower().endswith((".xls", ".xlsx")):
    M = pd.read_excel(args.matrix)
else:
    M = pd.read_csv(args.matrix, sep=None, engine="python")

first_col = M.columns[0]
M = M.rename(columns={first_col: "gene"}).set_index("gene")

# filtra genes no informativos
M_var = M[(M.sum(axis=1) > 0) & (M.sum(axis=1) < M.shape[1])]
if M_var.empty:
    raise SystemExit("[INFO] No hay variación (todo 0 o todo 1).")

# prioriza genes TNBC
tnbc_core = ["TP53","BRCA1","BRCA2","PIK3CA","PTEN","RB1","ARID1A"]
tnbc_in_matrix = [g for g in tnbc_core if g in M_var.index]

row_std = M_var.std(axis=1)
top = row_std.sort_values(ascending=False).index.tolist()

genes = []
if args.tnbc_first and tnbc_in_matrix:
    genes.extend(tnbc_in_matrix)
for g in top:
    if g not in genes:
        genes.append(g)
    if len(genes) >= args.topk:
        break

Mplot = M_var.loc[genes]

# heatmap
sns.set_context("talk")
plt.figure(figsize=(1.1*Mplot.shape[1], 0.24*Mplot.shape[0]))
ax = sns.heatmap(Mplot, cmap="Reds", cbar=False, linewidths=0.4, linecolor="lightgray")
ax.set_xlabel("Muestras"); ax.set_ylabel("Genes")
plt.title("Mutaciones somáticas (0/1) – genes más variables")
plt.tight_layout()
plt.savefig(args.out_png, dpi=300)
print(f"[OK] {args.out_png}")
