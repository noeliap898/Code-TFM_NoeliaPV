
import re, sys
import pandas as pd

def norm_sample(name: str) -> str:
    """
    Normaliza IDs de muestra quitando cualquier sufijo de versión que empiece por '~v'.
    - '...-organoid~v2'                 -> '...-organoid'
    - '...-organoid~v2.0.2.51.0'        -> '...-organoid'
    Además, recorta espacios.
    """
    name = str(name).strip()
    # quita todo desde ~v... hasta el final
    name = re.sub(r"~v[0-9].*$", "", name)
    return name

def load_tsv_index_rows(path):
    return pd.read_csv(path, sep="\t", index_col=0)

def harmonize_columns(df):
    df = df.copy()
    df.columns = [norm_sample(c) for c in df.columns]
    return df

def main(tpm_path, mut_path, out_tpm, out_mut):
    tpm = load_tsv_index_rows(tpm_path)
    mut = load_tsv_index_rows(mut_path)

    # normaliza índices (genes)
    tpm.index = tpm.index.astype(str).str.strip()
    mut.index = mut.index.astype(str).str.strip()

    # armoniza columnas (muestras)
    tpm = harmonize_columns(tpm)
    mut = harmonize_columns(mut)

    # busca intersección exacta
    common = [c for c in tpm.columns if c in mut.columns]
    if not common:
        print("[ERROR] No hay muestras en común tras armonizar. Revisa las columnas:")
        print("Ejemplo cols TPM:", list(tpm.columns)[:5])
        print("Ejemplo cols MUT:", list(mut.columns)[:5])
        sys.exit(1)

    # recorta a la intersección, preservando el orden de TPM
    tpm = tpm[common]
    mut = mut[common]

    # guarda
    tpm.to_csv(out_tpm, sep="\t")
    mut.to_csv(out_mut, sep="\t")
    print(f"[OK] Guardado TPM armonizado -> {out_tpm}  (n muestras={len(common)})")
    print(f"[OK] Guardado MUT armonizado -> {out_mut}  (n muestras={len(common)})")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Uso:")
        print("  python harmonize_ids_and_save.py <tpm_tsv> <mut_tsv> <out_tpm> <out_mut>")
        sys.exit(1)
    tpm_path, mut_path, out_tpm, out_mut = sys.argv[1:5]
    main(tpm_path, mut_path, out_tpm, out_mut)
