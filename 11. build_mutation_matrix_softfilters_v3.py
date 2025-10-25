
import os, sys, re, argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Panel TNBC que queremos priorizar al principio del heatmap
TNBC_PANEL = ["TP53","BRCA1","BRCA2","PTEN","RB1","PIK3CA","NF1","ARID1A","KMT2D","PTPN11"]

# Columnas candidatas
GENE_CANDIDATES   = ["Hugo_Symbol","Gene","Gene_Symbol","Gene_Name","gene_name","SYMBOL"]
SAMPLE_CANDIDATES = ["Tumor_Sample_Barcode","Tumor_Sample","Sample","Sample_ID","Tumor_Sample"]

def find_first_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def looks_like_ensembl_gene(x: str) -> bool:
    return bool(re.match(r"^ENSG\d+", str(x)))

def derive_sample_from_filename(fn):
    base = os.path.basename(fn)
    if base.startswith("._"):  # ficheros fantasma macOS
        base = base[2:]
    # quita sufijo .maf(.gz)
    base = re.sub(r"\.maf(\.gz)?$", "", base, flags=re.IGNORECASE)
    # si contiene ~WES o ~RNA, nos quedamos con el prefijo largo (ID del organoide)
    # p.ej. "388244~064-R~V1-organoid~v2.0.2.51.0"
    m = re.search(r"^(.*?)(~WES|~RNA|\.somatic|\.)", base)
    return (m.group(1) if m else base)

def read_maf_safe(path):
    # Ignora líneas de metadatos que comienzan por '##'
    try:
        return pd.read_csv(path, sep="\t", dtype=str, comment="#", engine="c", on_bad_lines="skip")
    except Exception:
        return pd.read_csv(path, sep="\t", dtype=str, comment="#", engine="python", on_bad_lines="skip")

def normalize_genes(df, gene_col):
    """Devuelve una Serie con símbolos de gen limpios.
       - Si Hugo_Symbol es ENSG*, intenta mapear a HGNC_Approved_name.
       - Quita vacíos y UNKNOWN.
    """
    genes = df[gene_col].fillna("").astype(str).str.strip()

    # Si hay columna HGNC_Approved_name, úsala para sustituir ENSG* por símbolo
    if "HGNC_Approved_name" in df.columns:
        mask_ensg = genes.apply(looks_like_ensembl_gene)
        if mask_ensg.any():
            # reemplaza por HGNC_Approved_name cuando esté disponible
            rep = df.loc[mask_ensg, "HGNC_Approved_name"].fillna("").astype(str).str.strip()
            # sólo sustituir si hay algo; si está vacío, mantenemos el valor original
            genes.loc[mask_ensg & rep.ne("")] = rep.loc[mask_ensg & rep.ne("")]

    # elimina vacíos y marcadores inútiles
    genes = genes.replace({"": np.nan, "__UNKNOWN__": np.nan, "UNKNOWN": np.nan}).dropna()

    # quita cualquier espacio y normaliza
    genes = genes.str.replace(r"\s+", "", regex=True)

    # por seguridad, eliminamos los que siguen siendo ENSG*
    genes = genes[~genes.apply(looks_like_ensembl_gene)]

    return genes

def collect_gene_sample_pairs(maf_dir):
    pairs = []  # (gene, sample)
    files = sorted([p for p in os.listdir(maf_dir)
                    if p.endswith(".maf") and not p.startswith("._")])
    if not files:
        print(f"[ERROR] No se encontraron .maf en {maf_dir}")
        sys.exit(1)

    for fname in files:
        fpath = os.path.join(maf_dir, fname)
        try:
            df = read_maf_safe(fpath)
        except Exception as e:
            print(f"[WARN] No pude leer {fname}: {e}")
            continue

        gene_col   = find_first_col(df, GENE_CANDIDATES)
        sample_col = find_first_col(df, SAMPLE_CANDIDATES)

        # Deriva SIEMPRE el sample_id del nombre del archivo, salvo que el MAF traiga algo válido distinto de __UNKNOWN__
        file_sample = derive_sample_from_filename(fname)
        if sample_col and df[sample_col].notna().any():
            any_valid = df[sample_col].fillna("").replace("__UNKNOWN__", "").str.strip()
            sample_id = any_valid.iloc[0] if (any_valid != "").any() else file_sample
        else:
            sample_id = file_sample

        # si no hay columna de gen, marcamos UNKNOWN (se filtrará)
        if gene_col is None:
            df["__gene__"] = "UNKNOWN"
            gene_col = "__gene__"

        # soft filters: prioriza somáticas si hay columna Mutation_Status
        df2 = df.copy()
        if "Mutation_Status" in df.columns:
            som = df["Mutation_Status"].fillna("").str.upper().str.contains("SOMATIC")
            if som.any():
                df2 = df.loc[som].copy()

        genes = normalize_genes(df2, gene_col)
        unique_genes = sorted(set(genes.tolist()))
        for g in unique_genes:
            pairs.append((g, sample_id))

    return pairs

def make_matrix(pairs):
    if not pairs:
        print("[ERROR] No se obtuvo ningún par gen×muestra (¿MAFs vacíos o todo ENSG/UNKNOWN?).")
        sys.exit(1)
    genes   = sorted(set([g for g,_ in pairs]))
    samples = sorted(set([s for _,s in pairs]))
    M = pd.DataFrame(0, index=genes, columns=samples, dtype=int)
    for g,s in pairs:
        M.loc[g, s] = 1
    return M

def reorder_with_tnbc_panel(M, top_n):
    # top por frecuencia
    freq = M.sum(axis=1).sort_values(ascending=False)
    top = list(freq.index[:top_n])

    # Panel TNBC (al principio si existe)
    tnbc_present = [g for g in TNBC_PANEL if g in M.index]
    ordered = tnbc_present + [g for g in top if g not in tnbc_present]
    ordered = list(dict.fromkeys(ordered))  # quita duplicados conservando orden
    if not ordered:
        return M
    return M.loc[ordered]

def plot_heatmap(M, out_png):
    if M.empty:
        print("[WARN] Matriz vacía; no hay figura.")
        return
    plt.figure(figsize=(max(8, 0.6*len(M.columns)), max(6, 0.35*len(M.index))))
    sns.heatmap(M, cmap="Reds", cbar=False, linewidths=0.4, linecolor="lightgray")
    plt.xlabel("Muestras (organoides)")
    plt.ylabel("Genes")
    plt.title("Mutaciones somáticas (presencia/ausencia)\nPanel TNBC priorizado + genes más frecuentes")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"[OK] Figura -> {out_png}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("maf_dir", help="Directorio con los .maf de Funcotator (uno por muestra)")
    ap.add_argument("out_dir", help="Directorio de salida")
    ap.add_argument("top_n", type=int, help="Número total de genes a mostrar (top_n)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Limpia ficheros '._*' (por si acaso)
    for junk in [p for p in os.listdir(args.maf_dir) if p.startswith("._")]:
        try:
            os.remove(os.path.join(args.maf_dir, junk))
        except Exception:
            pass

    pairs = collect_gene_sample_pairs(args.maf_dir)
    M = make_matrix(pairs)

    full_tsv = os.path.join(args.out_dir, "mutated_genes_matrix_full.tsv")
    M.to_csv(full_tsv, sep="\t")
    print(f"[OK] Matriz completa -> {full_tsv}")

    Msub = reorder_with_tnbc_panel(M, args.top_n)
    top_tsv = os.path.join(args.out_dir, "mutated_genes_matrix_top.tsv")
    Msub.to_csv(top_tsv, sep="\t")
    print(f"[OK] Matriz top -> {top_tsv}")

    out_png = os.path.join(args.out_dir, "figura3_heatmap_mutaciones.png")
    plot_heatmap(Msub, out_png)

if __name__ == "__main__":
    main()
