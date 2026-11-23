"""
MDR gap vs MIC shift across species (the paper's "Relationship Between MDR Differences and MIC Shifts").

For each species with MDR data, and each comparison with children (0-18 vs 19-64, 0-18 vs 65+):
  x = country-adjusted MDR difference (older minus 0-18, percentage points)
  y = median log2 MIC shift across every drug tested for that species and comparison (comparisons with >= 3
      qualifying countries, robust or not)
Species need >= 3 such drugs. Spearman rank correlation across species, separately per comparison.

Reads the outputs of the two pipelines (run them first); their locations are set in correlation_input.txt:
  padi_all_comparisons.csv   (padi_full_pipeline.py)
  mdr_all_comparisons.csv    (mdr_full_pipeline.py)
Writes, in results_folder:
  mdr_mic_species.csv       one row per species and comparison
  mdr_mic_correlation.csv   Spearman rho, p and number of species per comparison

Usage:
  python3 mdr_mic_correlation.py             (reads correlation_input.txt)
  python3 mdr_mic_correlation.py my_input.txt
"""
import argparse
import os

import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
COMPARISONS = ["0-18_vs_19-64", "0-18_vs_65+"]
MIN_DRUGS = 3


def run(mic_table, mdr_table, out_dir):
    mic = pd.read_csv(mic_table)
    mdr = pd.read_csv(mdr_table)
    rows = []
    for comp in COMPARISONS:
        m = mdr[(mdr.comparison == comp) & (mdr.phenotype_flag == "MDR") & (mdr.status != "insufficient_countries")]
        m = m.dropna(subset=["country_adjusted_rate_diff_pp"])
        for _, x in m.iterrows():
            d = mic[(mic.species == x.species) & (mic.comparison == comp) & (mic.status != "insufficient_countries")]
            if len(d) >= MIN_DRUGS:
                rows.append({"comparison": comp, "species": x.species, "mdr_gap_pp": x.country_adjusted_rate_diff_pp,
                             "mdr_robust": x.status == "robust", "median_log2_shift": d.log2_shift.median(),
                             "n_drugs": len(d)})
    df = pd.DataFrame(rows)
    cor = []
    for comp in COMPARISONS:
        d = df[df.comparison == comp]
        rho, p = spearmanr(d.mdr_gap_pp, d.median_log2_shift)
        cor.append({"comparison": comp, "rho": round(float(rho), 3), "p_value": float(p), "n_species": len(d)})
    cor = pd.DataFrame(cor)
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(os.path.join(out_dir, "mdr_mic_species.csv"), index=False)
    cor.to_csv(os.path.join(out_dir, "mdr_mic_correlation.csv"), index=False)
    print(cor.to_string(index=False))
    print(f"Wrote mdr_mic_species.csv and mdr_mic_correlation.csv to {out_dir}")
    return cor


# ============================================================================
# Input file. Each script has its own (see the instructions at the top of each one); the same reader is
# copied into all five scripts so each one runs on its own.
# ============================================================================
def norm_drug(d):
    return "".join(ch for ch in d.lower() if ch.isalnum())


def read_settings(path, need_bacteria=False):
    """Reads an input file: 'key = value' lines, plus (for the analysis steps) one line per bacterium to analyze,
    either 'Bacteria, Drug', 'Bacteria, all', 'Bacteria' alone (= every drug) or 'all' / 'all, all' (everything).
    Returns a dict of settings (paths resolved relative to the file) plus "pairs":
    {bacterium: set of normalized drug names, or None for every drug}, or None for every bacterium and drug."""
    if not os.path.exists(path):
        raise SystemExit(f"Settings file not found: {path}")
    base = os.path.dirname(os.path.abspath(path))
    raw, pairs, everything = {}, {}, False
    for n, line in enumerate(open(path), 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" in line:
            k, v = (x.strip() for x in line.split("=", 1))
            raw[k.lower()] = v
            continue
        sp, dr = (x.strip() for x in line.split(",", 1)) if "," in line else (line, "all")
        if sp.lower() == "all":
            everything = True
        elif dr.lower() == "all":
            pairs[sp] = None
        elif pairs.get(sp, set()) is not None:
            pairs.setdefault(sp, set()).add(norm_drug(dr))

    def folder(v):
        return v if os.path.isabs(v) else os.path.join(base, v)
    s = {"path": os.path.abspath(path)}
    for key in ("pediatric_folder", "age_group_folder", "mdr_dtr_folder", "year_folder", "results_folder",
                "sample_size_table", "padi_age_trends", "padi_all_comparisons", "mdr_all_comparisons"):
        s[key] = folder(raw[key]) if raw.get(key) else None
    s["specimen_source_folders"] = [folder(x.strip()) for x in raw.get("specimen_source_folders", "").split(",") if x.strip()]
    s["results_folder"] = s["results_folder"] or folder("results")
    s["permutations"] = int(raw.get("permutations", 0) or 0)
    s["workers"] = int(raw.get("workers", 4) or 4)
    s["pairs"] = None if everything else pairs
    if need_bacteria and not everything and not pairs:
        raise SystemExit(f"{path}: no bacteria listed (add a line such as 'all')")
    return s


def require_folder(s, key, step):
    """The folder setting `key`, or a clear error if it is missing."""
    d = s.get(key)
    if not d or not os.path.isdir(d):
        raise SystemExit(f"{key} not found: {d}\nSet {key} in {s['path']} to the folder holding the {step} files "
                         f"(the instructions at the top of that file show how to download and name them).")
    return d


def require_file(s, key, made_by):
    """The file setting `key`, or a clear error if it is missing."""
    p = s.get(key)
    if not p or not os.path.exists(p):
        raise SystemExit(f"{key} not found: {p}\nRun {made_by} first, or set {key} in {s['path']} to where its "
                         f"output is.")
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", default=os.path.join(HERE, "correlation_input.txt"), help="input file (default: correlation_input.txt)")
    args = ap.parse_args()
    s = read_settings(args.input)
    run(require_file(s, "padi_all_comparisons", "padi_full_pipeline.py"),
        require_file(s, "mdr_all_comparisons", "mdr_full_pipeline.py"), s["results_folder"])


if __name__ == "__main__":
    main()
