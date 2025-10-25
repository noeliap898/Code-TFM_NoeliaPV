

import os, re, glob, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def is_hidden_or_macdot(path: str) -> bool:
    base = os.path.basename(path)
    return base.startswith("._") or base.startswith(".")

def parse_stats_any(path: str) -> dict:
    """
    Devuelve un diccionario {filtro: conteo} a partir de un archivo
    *.filteringStats.tsv de GATK FilterMutectCalls en distintos formatos.
    """
    #Intento directo TSV/tab con columnas
    try:
        df = pd.read_csv(path, sep="\t")
        # normalización de nombres
        cols = [c.strip().lower() for c in df.columns]
        df.columns = cols
        # opciones:
        # a) columnas: filter | count
        if "filter" in df.columns and "count" in df.columns:
            d = {}
            for _, r in df.iterrows():
                key = str(r["filter"]).strip()
                try:
                    val = int(str(r["count"]).strip())
                except Exception:
                    val = pd.to_numeric(r["count"], errors="coerce")
                    val = 0 if pd.isna(val) else int(val)
                if key:
                    d[key] = d.get(key, 0) + val
            if d:
                return d
        # b) dos primeras columnas ya son filtro/numero sin cabecera clara
        if df.shape[1] >= 2:
            d = {}
            for _, r in df.iterrows():
                key = str(r.iloc[0]).strip()
                val = pd.to_numeric(r.iloc[1], errors="coerce")
                if key and not pd.isna(val):
                    d[key] = d.get(key, 0) + int(val)
            if d:
                return d
    except Exception:
        pass

    # parseo “línea a línea”: soporta "key: value", "key = value",
    #    "key\tvalue" e ignora cabeceras.
    d = {}
    with open(path, "r", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # elimina cabeceras tipo "filter count", "Metric  Count", etc.
            if re.match(r"(?i)^(filter|metric|name)\s+(\t|,|\s)+count", line):
                continue

            # separadores posibles
            if "\t" in line:
                parts = [p.strip() for p in line.split("\t") if p.strip() != ""]
            elif ":" in line:
                parts = [p.strip() for p in line.split(":", 1)]
            elif "=" in line:
                parts = [p.strip() for p in line.split("=", 1)]
            else:
                # líneas que no encajan (las ignoramos)
                continue

            if len(parts) < 2:
                continue

            key, val_raw = parts[0], parts[1]
            val_raw = re.split(r"[^\d]+", val_raw.strip())[0] if re.search(r"\d", val_raw) else ""
            if not val_raw:
                continue

            try:
                val = int(val_raw)
            except Exception:
                val = pd.to_numeric(val_raw, errors="coerce")
                if pd.isna(val):
                    continue
                val = int(val)

            if key:
                d[key] = d.get(key, 0) + val

    return d

def clean_sample_name(fname: str) -> str:
    # elimina prefijo "._" si existe
    base = os.path.basename(fname).lstrip("._")
    # quita el sufijo a partir de "~WES.somatic.filtered.vcf.gz.filteringStats.tsv"
    base = re.sub(r"~WES\.somatic.*?\.filteringStats\.tsv$", "", base)
    return base

def build_metrics(in_dir: str) -> pd.DataFrame:
    paths = sorted(glob.glob(os.path.join(in_dir, "*.filteringStats.tsv")))
    rows = []
    for p in paths:
        if is_hidden_or_macdot(p):
            continue
        d = parse_stats_any(p)
        if not d:
            # si no podemos parsear, lo anotamos, pero seguimos
            print(f"[WARN] No se pudieron extraer conteos de: {p}")
            continue
        sample = clean_sample_name(p)

        # normaliza claves (p. ej., "PASS"/ "Pass")
        nd = {}
        for k, v in d.items():
            nd[str(k).strip().upper()] = nd.get(str(k).strip().upper(), 0) + int(v)

        pass_count = nd.get("PASS", 0)
        filtered_total = sum(v for k, v in nd.items() if k != "PASS")

        row = {"sample": sample, "PASS": pass_count, "Filtered_Total": filtered_total}
        # añade motivos individuales (columnas) para exploración
        for k, v in nd.items():
            if k == "PASS":
                continue
            row[f"FILTER_{k}"] = v
        rows.append(row)

    if not rows:
        return pd.DataFrame(columns=["sample", "PASS", "Filtered_Total"])

    df = pd.DataFrame(rows)
    # orden por PASS descendente (más informativo)
    df = df.sort_values(["PASS", "Filtered_Total"], ascending=[False, True]).reset_index(drop=True)
    return df

def plot_stacked_pass_filtered(df: pd.DataFrame, out_png: str):
    # reemplaza NaN por 0
    P = df.copy()
    for c in ["PASS", "Filtered_Total"]:
        if c not in P.columns: P[c] = 0
    P[["PASS", "Filtered_Total"]] = P[["PASS", "Filtered_Total"]].fillna(0).astype(int)

    x = np.arange(len(P))
    plt.figure(figsize=(12, 5))
    plt.bar(x, P["PASS"], label="PASS")
    plt.bar(x, P["Filtered_Total"], bottom=P["PASS"], label="Filtradas")

    plt.xticks(x, P["sample"], rotation=45, ha="right")
    plt.ylabel("Nº variantes")
    plt.title("Mutect2 FilterMutectCalls: PASS vs Filtradas")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_percent_pass(df: pd.DataFrame, out_png: str):
    P = df.copy()
    for c in ["PASS", "Filtered_Total"]:
        if c not in P.columns: P[c] = 0
    P[["PASS", "Filtered_Total"]] = P[["PASS", "Filtered_Total"]].fillna(0).astype(int)
    P["Total"] = P["PASS"] + P["Filtered_Total"]
    # evita división por cero
    P["Percent_PASS"] = np.where(P["Total"] > 0, 100.0 * P["PASS"] / P["Total"], 0.0)

    plt.figure(figsize=(12, 3.8))
    plt.bar(np.arange(len(P)), P["Percent_PASS"])
    plt.xticks(np.arange(len(P)), P["sample"], rotation=45, ha="right")
    plt.ylabel("% PASS")
    plt.ylim(0, 100)
    plt.title("Porcentaje de variantes que pasan filtro")
    # si todo es 0, lo indicamos en el subtítulo
    if (P["Total"] == 0).all():
        plt.suptitle("Aviso: todos los totales son 0 (¿ficheros vacíos?)", y=0.97, fontsize=9)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def main():
    if len(sys.argv) < 3:
        print("Uso:")
        print("  python wes_filtering_qc_v2.py <dir_filtrados> <dir_salida>")
        print("\nEjemplo:")
        print("  python wes_filtering_qc_v2.py /Volumes/SecTFM/results/wes /Volumes/SecTFM/results/figuras")
        sys.exit(0)

    in_dir  = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    df = build_metrics(in_dir)
    out_tsv = os.path.join(out_dir, "wes_filtering_metrics.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[OK] Guardado: {out_tsv}")

    plot_stacked_pass_filtered(df, os.path.join(out_dir, "fig2a_wes_pass_vs_filtered_stacked.png"))
    plot_percent_pass(df,        os.path.join(out_dir, "fig2b_wes_percent_pass.png"))
    print(f"[OK] Figuras en: {out_dir}")

if __name__ == "__main__":
    main()
