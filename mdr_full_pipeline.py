"""
MDR/DTR full pipeline: every step of the paper's MDR/DTR age analysis in one run.

Give it species. Assuming the ATLAS exports are already on disk, it outputs how the share of multidrug-resistant
(MDR) and difficult-to-treat resistant (DTR) isolates differs between age groups, and which differences are robust.

Steps, in order (same math and cutoffs as the paper):
  1. For each species and age group, read three exports (View By: Country, all drugs): no phenotype filter
     (baseline), MDR (CLSI) filter, DTR filter. Isolate counts come from each export's "Species" section.
  2. Per country: %MDR = N(MDR) / N(baseline), %DTR = N(DTR) / N(baseline). A country qualifies when its baseline
     has >= 20 isolates in both age groups; a country missing from the MDR/DTR export counts as 0 (ATLAS leaves
     out countries with no matching isolates). >= 3 qualifying countries needed.
  3. Country-adjusted difference (older minus younger, percentage points), each country weighted by the harmonic
     mean of its two baseline counts. Country bootstrap (20,000 draws, binomial resampling, fixed seed) -> 95% CI
     and two-sided p-value.
  4. Benjamini-Hochberg across every MDR/DTR comparison tested in this run (one family, separate from the MIC
     shifts). Robust = q < 0.05. No effect size floor.
  5. Testing panel check: per country and age group, coverage of each drug class = the highest share of isolates
     tested against any one drug in that class; the coverage index is the sum over classes. Panels are
     "comparable" when the indices differ by < 0.5 classes and no class differs by >= 10 percentage points.
     Reported, not gating.
  6. Specimen source and year composition of the baseline exports compared between age groups; flagged when
     TVD > 15 pp and chi-square p < 0.01.
  7. For robust, flagged comparisons: the same difference recomputed within each specimen source (or year) from
     worldwide single-source (single-year) exports, weighted and bootstrapped (2,000 draws). The finding survives
     if that CI excludes zero; with no usable source (year) it is "not checked". robust_final = survives (or did
     not need) both checks.

The BH family is every comparison tested in THIS run, so q-values match the paper only when the run covers the
same species (the default: every species folder in the data directory).

Inputs: set in mdr_input.txt (mdr_dtr_folder, specimen_source_folders, year_folder, results_folder), with the
bacteria to analyze ("all" = every bacterium in mdr_dtr_folder).
Outputs, in results_folder:
  mdr_all_comparisons.csv   every species x flag x comparison computed, with its status
  mdr_trends.csv            robust comparisons with rates, difference, CI, q, panel check, specimen source and year
  mdr_summary.json          counts at each step
  files_to_download.txt     only if a flagged finding still needs specimen source or year files

Usage:
  python3 mdr_full_pipeline.py                  (reads mdr_input.txt)
  python3 mdr_full_pipeline.py my_input.txt
"""
import argparse
import csv
import json
import os
import re

import numpy as np
import openpyxl
from scipy.stats import chi2_contingency

# ============================================================================
# Settings (the paper's values)
# ============================================================================
MIN_COUNTRY_N = 20                # baseline isolates per country (or source / year) per age group
MIN_QUALIFYING_COUNTRIES = 3
N_BOOTSTRAP = 2000                # specimen source / year checks and the panel check
N_BOOTSTRAP_GATE = 20000          # main country bootstrap
FDR_Q = 0.05
RNG_SEED = 0
AGE_GROUPS = ["0-18", "19-64", "65+"]
COMPARISONS = [("0-18", "19-64"), ("19-64", "65+"), ("0-18", "65+")]
FLAGS = [("MDR", "mdr"), ("DTR", "dtr")]
INDEX_DIFF_LIMIT = 0.5            # testing panel: coverage index difference, classes
CLASS_DIFF_LIMIT = 10.0           # testing panel: one class's coverage difference, percentage points
TITLE_RE = re.compile(r"from (.+?)(?:\s*\(N|$)")
PANEL_NAME = {"Enterobacter hormaechei": "Enterobacter hormaechi"}   # how ATLAS spells it in the Species section

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA_DIR = os.path.join(HERE, "data", "mdr_dtr")
DEFAULT_SITE_DIR = os.path.join(HERE, "data", "specimen_source")
DEFAULT_SITE_MANIFESTS = []
DEFAULT_YEAR_DIR = os.path.join(HERE, "data", "year")
DEFAULT_OUT_DIR = os.path.join(HERE, "results", "mdr_dtr")

# Drug -> class, keyed by lower-case letters and digits only (so 'Pip/taz' and 'Piperacillin-tazobactam' both match)
DRUG_CLASS = {
    'amikacin': 'Aminoglycoside',
    'amoxicillinclavulanate': 'Penicillin',
    'amoxyclav': 'Penicillin',
    'ampicillin': 'Penicillin',
    'ampicillinsulbactam': 'Penicillin',
    'azithromycin': 'Macrolide',
    'aztreonam': 'Monobactam',
    'aztreonamavibactam': 'Monobactam',
    'cefepime': 'Cephalosporin',
    'cefiderocol': 'Cephalosporin',
    'cefixime': 'Cephalosporin',
    'cefoperazonesulbactam': 'Cephalosporin',
    'cefoxitin': 'Cephalosporin',
    'cefpodoxime': 'Cephalosporin',
    'ceftaroline': 'Cephalosporin',
    'ceftarolineavibactam': 'Cephalosporin',
    'ceftazidime': 'Cephalosporin',
    'ceftazidimeavibactam': 'Cephalosporin',
    'ceftibuten': 'Cephalosporin',
    'ceftibutenavibactam': 'Cephalosporin',
    'ceftolozanetazobactam': 'Cephalosporin',
    'ceftriaxone': 'Cephalosporin',
    'ciprofloxacin': 'Quinolone',
    'clarithromycin': 'Macrolide',
    'clindamycin': 'Lincosamide',
    'colistin': 'Miscellaneous',
    'daptomycin': 'Lipopeptides',
    'doripenem': 'Carbapenem',
    'ertapenem': 'Carbapenem',
    'erythromycin': 'Macrolide',
    'gatifloxacin': 'Quinolone',
    'gentamicin': 'Aminoglycoside',
    'imipenem': 'Carbapenem',
    'levofloxacin': 'Quinolone',
    'linezolid': 'Miscellaneous',
    'meropenem': 'Carbapenem',
    'meropenemvaborbactam': 'Carbapenem',
    'metronidazole': 'Miscellaneous',
    'minocycline': 'Tetracycline',
    'moxifloxacin': 'Quinolone',
    'oxacillin': 'Penicillin',
    'penicillin': 'Penicillin',
    'piperacillintazobactam': 'Penicillin',
    'piptaz': 'Penicillin',
    'quinupristindalfopristin': 'Streptogramins',
    'sulbactam': 'Miscellaneous',
    'tebipenem': 'Carbapenem',
    'teicoplanin': 'Glycopeptide',
    'tetracycline': 'Tetracycline',
    'tigecycline': 'Glycylcycline',
    'tobramycin': 'Aminoglycoside',
    'trimethoprimsulfa': 'Folate inhibitor',
    'trimethoprimsulfamethoxazole': 'Folate inhibitor',
    'vancomycin': 'Glycopeptide',
}


# ---- parsing and shared math (the same functions as in padi_full_pipeline.py)
def parse_by_country_blocks(path, sheet="Absolute"):
    """{country: {drug: {"N": int, "counts": {mic: count}}}}"""
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(min_row=1, max_row=ws.max_row, values_only=True))
    data, current, header, in_table = {}, None, None, False
    for row in rows:
        first = row[0]
        if first and isinstance(first, str) and ("isolates at each MIC" in first or "frequency distribution" in first):
            m = TITLE_RE.search(first)
            current = m.group(1).strip() if m else first
            data.setdefault(current, {})
            header, in_table = None, False
            continue
        if first == "Antibacterial":
            header, in_table = row, True
            continue
        if not in_table or current is None:
            continue
        if first is None:
            in_table = False
            continue
        if not isinstance(first, str) or not isinstance(row[1], (int, float)):
            in_table = False
            continue
        counts = {}
        for mic_val, c in zip(header[2:], row[2:]):
            if mic_val is None or c is None:
                continue
            counts[float(mic_val)] = int(c)
        data[current][first] = {"N": int(row[1]), "counts": counts}
    return data


def parse_composition_sections(path, sheet="Absolute"):
    """{country: {section_name: {category: frequency}}} -- Specimen Source
    Group / Ward / Year / Phenotype metadata already embedded ahead of the
    MIC table in every by-country export, free of charge."""
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(min_row=1, max_row=ws.max_row, values_only=True))
    data, current_country, current_section = {}, None, None
    i, n = 0, len(rows)
    while i < n:
        row = rows[i]
        first = row[0]
        if first and isinstance(first, str) and "frequency distribution" in first:
            m = TITLE_RE.search(first)
            current_country = m.group(1).strip() if m else first
            data.setdefault(current_country, {})
            current_section = None
        elif (isinstance(first, str) and row[1] is None and row[2] is None
              and i + 1 < n and rows[i + 1][:3] == (first, "Frequency", "Percentage of Total")):
            current_section = first
            data[current_country].setdefault(current_section, {})
            i += 1
        elif current_section and current_country and first is not None and isinstance(row[1], (int, float)):
            data[current_country][current_section][first] = data[current_country][current_section].get(first, 0) + int(row[1])
        i += 1
    return data


def pooled_category_counts(by_country_composition, section):
    pooled = {}
    for sections in by_country_composition.values():
        for cat, freq in sections.get(section, {}).items():
            pooled[cat] = pooled.get(cat, 0) + freq
    return pooled


def compare_composition(counts_a, counts_b):
    cats = sorted(set(counts_a) | set(counts_b))
    total_a, total_b = sum(counts_a.values()), sum(counts_b.values())
    if total_a == 0 or total_b == 0:
        return None
    tvd = 0.5 * sum(abs(counts_a.get(c, 0) / total_a - counts_b.get(c, 0) / total_b) for c in cats)
    table = np.array([[counts_a.get(c, 0) for c in cats], [counts_b.get(c, 0) for c in cats]])
    table = table[:, table.sum(axis=0) > 0]
    if table.shape[1] < 2:
        return {"tvd_pp": round(tvd * 100, 1), "chi2_p": None}
    chi2, p, dof, _ = chi2_contingency(table)
    return {"tvd_pp": round(tvd * 100, 1), "chi2_p": float(p)}


def bh_qvalues(pvals):
    """Benjamini-Hochberg adjusted p-values (q-values), same order as input. None stays None and is NOT
    counted in m -- callers must pass every test in the family."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    q = [None] * len(pvals)
    m = len(idx)
    if m == 0:
        return q
    p = np.array([pvals[i] for i in idx], dtype=float)
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.clip(ranked, 0, 1)
    for r, o in enumerate(order):
        q[idx[o]] = float(ranked[r])
    return q



# ---- rate difference and bootstrap
def weighted_rate_diff(entries):
    """entries: [(unit, n_a, pos_a, n_b, pos_b)], where a unit is a country, specimen source or year.
    Same harmonic-mean-N weighting as the MIC shift."""
    diffs, weights = [], []
    for _, n_a, pos_a, n_b, pos_b in entries:
        ra, rb = pos_a / n_a, pos_b / n_b
        diffs.append(rb - ra)
        weights.append(2 * n_a * n_b / (n_a + n_b))
    if not diffs:
        return None
    return sum(w * d for w, d in zip(weights, diffs)) / sum(weights)


def bootstrap_rate_ci(entries, n_boot=N_BOOTSTRAP, seed=RNG_SEED):
    """Cluster bootstrap over units (countries, specimen sources or years): resample which units, with
    replacement, then binomial-resample each unit's positive count given its fixed N and observed rate.
    Same structure as the MIC bootstrap, with a binomial instead of a multinomial because MDR/DTR (or a marker)
    is a yes/no outcome per isolate, not a MIC category."""
    rng = np.random.default_rng(seed)
    n = len(entries)
    if n == 0:
        return None
    # Vectorized 2026-11-08 (same design: resample units, then binomial per pick) + bootstrap p-value.
    n_arr = np.array([[e[1], e[3]] for e in entries], dtype=float)
    pos_arr = np.array([[e[2], e[4]] for e in entries], dtype=float)
    idx = rng.integers(0, n, size=(n_boot, n))
    na, nb = n_arr[idx, 0], n_arr[idx, 1]
    ra = rng.binomial(na.astype(int), pos_arr[idx, 0] / na) / na
    rb = rng.binomial(nb.astype(int), pos_arr[idx, 1] / nb) / nb
    w = 2 * na * nb / (na + nb)
    arr = (w * (rb - ra)).sum(axis=1) / w.sum(axis=1)
    lo, hi = np.percentile(arr, [2.5, 97.5])
    point = weighted_rate_diff(entries)
    far = int((arr <= 0).sum()) if point > 0 else int((arr >= 0).sum()) if point < 0 else arr.size
    return {"ci_lo_pp": round(float(lo) * 100, 2), "ci_hi_pp": round(float(hi) * 100, 2),
            "excludes_zero": bool(lo > 0 or hi < 0),
            "p_value": min(1.0, 2.0 * (far + 1) / (arr.size + 1)), "n_boot": int(n_boot)}



# ============================================================================
# Helpers
# ============================================================================
def safe(s):
    return s.replace(" ", "_").replace(",", "")


def safe_site(s):
    return s.replace(" ", "_").replace("/", "-")


def norm_name(d):
    return "".join(ch for ch in d.lower() if ch.isalnum())


def drug_class(drug):
    cls = DRUG_CLASS.get(norm_name(drug))
    if cls is None:
        raise SystemExit(f"no drug class for {drug!r}: add it to DRUG_CLASS")
    return cls


def species_n_by_country(path, species):
    """Isolate N per country from the export's 'Species' section. A country absent from the file is absent here;
    an explicit 0 stays 0."""
    out = {}
    for country, secs in parse_composition_sections(path).items():
        sp = secs.get("Species", {})
        if species in sp:
            out[country] = sp[species]
    return out


def species_n(path, species):
    """Isolate N from a worldwide (single specimen source or single year) export's 'Species' section, or None."""
    for _, secs in parse_composition_sections(path).items():
        sp = secs.get("Species", {})
        if species in sp:
            return sp[species]
    return None


def load_manifest(paths):
    m = {}
    for p in paths:
        if not os.path.exists(p):
            continue
        for sp, groups in json.load(open(p)).items():
            for g, levels in groups.items():
                have = m.setdefault(sp, {}).setdefault(g, [])
                have.extend(x for x in levels if x not in have)
    return m


# ============================================================================
# Steps
# ============================================================================
def country_comparison(baseline_n, filtered_n, a, b):
    """Country-adjusted MDR or DTR difference for one comparison, or None if no country qualifies."""
    countries = sorted(set(baseline_n[a]) & set(baseline_n[b]))
    qualifying = [c for c in countries if baseline_n[a][c] >= MIN_COUNTRY_N and baseline_n[b][c] >= MIN_COUNTRY_N]
    if not qualifying:
        return None
    pos_a = sum(filtered_n[a].get(c, 0) for c in qualifying)
    n_a = sum(baseline_n[a][c] for c in qualifying)
    pos_b = sum(filtered_n[b].get(c, 0) for c in qualifying)
    n_b = sum(baseline_n[b][c] for c in qualifying)
    if n_a == 0 or n_b == 0:
        return None
    table = np.array([[pos_a, n_a - pos_a], [pos_b, n_b - pos_b]])
    if (table < 0).any():
        return None
    _, p, _, _ = chi2_contingency(table)
    diffs, weights, ra_vals, rb_vals, n_up, n_down = [], [], [], [], 0, 0
    for c in qualifying:
        ra = filtered_n[a].get(c, 0) / baseline_n[a][c]
        rb = filtered_n[b].get(c, 0) / baseline_n[b][c]
        w = 2 * baseline_n[a][c] * baseline_n[b][c] / (baseline_n[a][c] + baseline_n[b][c])
        diffs.append(rb - ra); weights.append(w); ra_vals.append(ra); rb_vals.append(rb)
        n_up += rb > ra
        n_down += rb < ra
    sw = sum(weights)
    boot = bootstrap_rate_ci([(c, baseline_n[a][c], filtered_n[a].get(c, 0), baseline_n[b][c], filtered_n[b].get(c, 0))
                              for c in qualifying], n_boot=N_BOOTSTRAP_GATE)
    return {
        "n_countries": len(qualifying),
        "pooled_rate_a_pct": round(pos_a / n_a * 100, 2), "pooled_rate_b_pct": round(pos_b / n_b * 100, 2),
        "pooled_n_a": n_a, "pooled_n_b": n_b, "chi2_p": float(p),
        "country_adjusted_rate_a_pct": round(sum(w * r for w, r in zip(weights, ra_vals)) / sw * 100, 2),
        "country_adjusted_rate_b_pct": round(sum(w * r for w, r in zip(weights, rb_vals)) / sw * 100, 2),
        "country_adjusted_rate_diff_pp": round(sum(w * d for w, d in zip(weights, diffs)) / sw * 100, 2),
        "bootstrap_ci_lo_pp": boot["ci_lo_pp"] if boot else "", "bootstrap_ci_hi_pp": boot["ci_hi_pp"] if boot else "",
        "bootstrap_p_value": boot["p_value"] if boot else "",
        "countries_higher_in_b": int(n_up), "countries_higher_in_a": int(n_down),
    }


def panel_coverage(baseline_path_a, baseline_path_b, species, rate_gaps):
    """Testing panel comparison between two age groups (see step 5)."""
    def coverage(path):
        totals = species_n_by_country(path, species)
        blocks = parse_by_country_blocks(path)
        out = {}
        for c, n_tot in totals.items():
            if not n_tot:
                continue
            cov = {}
            for drug, e in (blocks.get(c) or {}).items():
                cls = drug_class(drug)
                cov[cls] = max(cov.get(cls, 0.0), min(1.0, e["N"] / n_tot))
            out[c] = (n_tot, cov)
        return out
    cov_a, cov_b = coverage(baseline_path_a), coverage(baseline_path_b)
    qual = [c for c in sorted(set(cov_a) & set(cov_b)) if cov_a[c][0] >= MIN_COUNTRY_N and cov_b[c][0] >= MIN_COUNTRY_N]
    if not qual:
        return {"panel_coverage_verdict": "not computed"}
    classes = sorted({k for c in qual for k in list(cov_a[c][1]) + list(cov_b[c][1])})
    w = np.array([2 * cov_a[c][0] * cov_b[c][0] / (cov_a[c][0] + cov_b[c][0]) for c in qual])
    ca = np.array([[cov_a[c][1].get(k, 0.0) for k in classes] for c in qual])
    cb = np.array([[cov_b[c][1].get(k, 0.0) for k in classes] for c in qual])
    idx_a, idx_b = ca.sum(axis=1), cb.sum(axis=1)
    idx_diff = float((w * (idx_b - idx_a)).sum() / w.sum())
    cls_diff = (w[:, None] * (cb - ca)).sum(axis=0) / w.sum() * 100
    rng = np.random.default_rng(RNG_SEED)
    pick = rng.integers(0, len(qual), size=(N_BOOTSTRAP, len(qual)))
    boot = (w[pick] * (idx_b - idx_a)[pick]).sum(axis=1) / w[pick].sum(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    worst = int(np.argmax(np.abs(cls_diff)))
    differs = abs(idx_diff) >= INDEX_DIFF_LIMIT or abs(cls_diff[worst]) >= CLASS_DIFF_LIMIT
    return {"panel_coverage_index_a": round(float((w * idx_a).sum() / w.sum()), 2),
            "panel_coverage_index_b": round(float((w * idx_b).sum() / w.sum()), 2),
            "panel_coverage_index_diff": round(idx_diff, 3),
            "panel_coverage_index_diff_ci": f"{lo:+.3f} to {hi:+.3f}",
            "panel_max_class_diff_pp": round(float(cls_diff[worst]), 1), "panel_max_class_diff_class": classes[worst],
            "panel_classes_differing_10pp": "; ".join(f"{k} {v:+.1f}" for k, v in zip(classes, cls_diff)
                                                     if abs(v) >= CLASS_DIFF_LIMIT),
            "panel_coverage_verdict": "PANEL DIFFERS" if differs else "comparable",
            "_idx_diff": idx_diff, "_differs": differs}


def panel_note(pc, a, b, flag, gap_pp):
    """Whether a testing gap could produce (or works against) this MDR/DTR difference."""
    if not pc.get("_differs"):
        return ""
    less = b if pc["_idx_diff"] < 0 else a
    lower_rate = a if gap_pp > 0 else b
    return (f"{flag}: {less} tested less and has the {'lower' if less == lower_rate else 'higher'} rate -> "
            f"{'could contribute to' if less == lower_rate else 'works against'} the gap")


def composition_flags(path_a, path_b):
    out = {}
    comp_a, comp_b = parse_composition_sections(path_a), parse_composition_sections(path_b)
    for section, key in [("Specimen Source Group", "specimen_source"), ("Year", "year")]:
        cmp = compare_composition(pooled_category_counts(comp_a, section), pooled_category_counts(comp_b, section))
        if cmp is None or cmp.get("chi2_p") is None:
            out[key] = {"tvd_pp": "", "differs": False}
        else:
            out[key] = {"tvd_pp": cmp["tvd_pp"], "differs": bool(cmp["tvd_pp"] > 15 and cmp["chi2_p"] < 0.01)}
    return out


def levels_from_files(folder, g):
    """Specimen sources (or years) that have an unfiltered export in folder: atlas_pair_<g>_<level>_baseline.xlsx."""
    if not os.path.isdir(folder):
        return []
    pre, suf = f"atlas_pair_{g}_", "_baseline.xlsx"
    return sorted(fn[len(pre):-len(suf)].replace("_-_", " / ").replace("_", " ")
                  for fn in os.listdir(folder) if fn.startswith(pre) and fn.endswith(suf))


def needed_files(section, kind, a, b, baseline_a, baseline_b, file_for):
    """Files a flagged comparison needs for its specimen source or year check: for every level with
    >= MIN_COUNTRY_N isolates in both age groups' unfiltered exports, the unfiltered and the MDR (or DTR) export.
    Returns [(age group, level, file kind, path, already there?)]."""
    ca = pooled_category_counts(parse_composition_sections(baseline_a), section)
    cb = pooled_category_counts(parse_composition_sections(baseline_b), section)
    levels = sorted(str(lv) for lv in set(ca) & set(cb) if ca[lv] >= MIN_COUNTRY_N and cb[lv] >= MIN_COUNTRY_N)
    out = []
    for lv in levels:
        for g in (a, b):
            for k in ("baseline", kind):
                p = file_for(g, lv, k)
                out.append((g, lv, k, p, os.path.exists(p)))
    return out


def stratified_rate_check(species, kind, a, b, levels_a, levels_b, file_for):
    """MDR or DTR difference recomputed within each level of a factor (specimen source or year).
    file_for(age_group, level, file_kind) -> path. A level is used when both age groups have its baseline and its
    MDR/DTR export, with >= MIN_COUNTRY_N baseline isolates in each."""
    def n_for(g, levels):
        out = {}
        for lv in levels:
            v = {}
            for k in ("baseline", kind):
                p = file_for(g, lv, k)
                if os.path.exists(p):
                    n = species_n(p, species)
                    if n is not None:
                        v[k] = n
            if "baseline" in v and kind in v:
                out[lv] = v
        return out
    la, lb = n_for(a, levels_a), n_for(b, levels_b)
    qual = [lv for lv in sorted(set(la) & set(lb)) if la[lv]["baseline"] >= MIN_COUNTRY_N and lb[lv]["baseline"] >= MIN_COUNTRY_N]
    if not qual:
        return {"check": "flagged, not checked (no usable exports)", "n_used": 0}, "unknown"
    entries = [(lv, la[lv]["baseline"], la[lv][kind], lb[lv]["baseline"], lb[lv][kind]) for lv in qual]
    diff = weighted_rate_diff(entries)
    boot = bootstrap_rate_ci(entries)
    survives = bool(boot and boot["excludes_zero"])
    per = "; ".join(f"{lv}: {(lb[lv][kind] / lb[lv]['baseline'] - la[lv][kind] / la[lv]['baseline']) * 100:+.1f}pp"
                    for lv in qual)
    return {"check": "survives adjustment" if survives else "weakened by adjustment (CI crosses zero)",
            "n_used": len(qual), "used": "; ".join(qual), "adjusted_diff_pp": round(diff * 100, 2),
            "ci_lo_pp": (boot or {}).get("ci_lo_pp"), "ci_hi_pp": (boot or {}).get("ci_hi_pp"),
            "per_level": per}, survives


# ============================================================================
# Main run
# ============================================================================
ALL_FIELDS = [
    "species", "phenotype_flag", "comparison", "status", "n_countries",
    "pooled_rate_a_pct", "pooled_rate_b_pct", "pooled_n_a", "pooled_n_b", "chi2_p",
    "country_adjusted_rate_a_pct", "country_adjusted_rate_b_pct", "country_adjusted_rate_diff_pp",
    "bootstrap_ci_lo_pp", "bootstrap_ci_hi_pp", "bootstrap_p_value", "bh_q_value",
    "countries_higher_in_b", "countries_higher_in_a",
    "panel_coverage_verdict", "panel_coverage_index_a", "panel_coverage_index_b", "panel_coverage_index_diff",
    "panel_coverage_index_diff_ci", "panel_max_class_diff_pp", "panel_max_class_diff_class",
    "panel_classes_differing_10pp", "panel_coverage_note",
    "specimen_source_tvd_pp", "specimen_source_differs", "year_tvd_pp", "year_differs",
    "specimen_source_check", "site_adjusted_rate_diff_pp", "site_adjusted_ci_lo_pp", "site_adjusted_ci_hi_pp",
    "n_sites_used", "sites_used", "per_site_diffs", "robust_after_specimen_source_check",
    "year_check", "year_adjusted_rate_diff_pp", "year_adjusted_ci_lo_pp", "year_adjusted_ci_hi_pp",
    "n_years_used", "years_used", "per_year_diffs", "robust_after_year_check", "robust_final",
]


def run(species_list=None, data_dir=DEFAULT_DATA_DIR, site_folders=None, site_manifests=(),
        year_dir=DEFAULT_YEAR_DIR, out_dir=DEFAULT_OUT_DIR):
    if species_list is None:
        species_list = sorted(d.replace("_", " ") for d in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, d)))
    print(f"MDR/DTR full pipeline: {len(species_list)} species", flush=True)
    rows, skipped, baseline_paths = [], [], {}
    for species in species_list:
        name = PANEL_NAME.get(species, species)
        d = os.path.join(data_dir, safe(species))
        paths = {(g, k): os.path.join(d, f"atlas_pair_{g}_{k}.xlsx") for g in AGE_GROUPS for k in ("baseline", "mdr", "dtr")}
        paths = {key: p for key, p in paths.items() if os.path.exists(p)}
        have = {kind: {g for g in AGE_GROUPS if (g, "baseline") in paths and (g, kind) in paths} for _, kind in FLAGS}
        if all(len(v) < 2 for v in have.values()):
            skipped.append(species)
            continue
        baseline_paths[species] = {g: paths[(g, "baseline")] for g in AGE_GROUPS if (g, "baseline") in paths}
        baseline_n = {g: species_n_by_country(p, name) for g, p in baseline_paths[species].items()}
        for flag, kind in FLAGS:
            if len(have[kind]) < 2:
                continue
            filtered_n = {g: species_n_by_country(paths[(g, kind)], name) for g in have[kind]}
            for a, b in COMPARISONS:
                if a not in have[kind] or b not in have[kind]:
                    continue
                res = country_comparison(baseline_n, filtered_n, a, b)
                if res is None:
                    continue
                row = {f: "" for f in ALL_FIELDS} | res | {"species": species, "phenotype_flag": flag,
                                                         "comparison": f"{a}_vs_{b}"}
                row["status"] = ("insufficient_countries" if res["n_countries"] < MIN_QUALIFYING_COUNTRIES else "tested")
                rows.append(row)
        print(f"  {species}: {sum(1 for r in rows if r['species'] == species)} comparisons computed", flush=True)

    # ---- Benjamini-Hochberg over every comparison with >= 3 countries
    rows.sort(key=lambda r: (r["species"], r["phenotype_flag"], r["comparison"]))
    fam = [r for r in rows if r["status"] == "tested" and r["bootstrap_p_value"] != ""]
    for r, q in zip(fam, bh_qvalues([float(r["bootstrap_p_value"]) for r in fam])):
        r["bh_q_value"] = q
        r["status"] = ("robust" if q < FDR_Q else
                       "fails_bh_fdr" if float(r["bootstrap_p_value"]) < FDR_Q else "not_significant")
    robust = [r for r in fam if r["status"] == "robust"]
    print(f"BH family m = {len(fam)}; robust (q < {FDR_Q}): {len(robust)}", flush=True)

    # ---- Panel check and composition flags (every computed comparison), then factor checks (robust + flagged)
    panel_cache, comp_cache = {}, {}
    manifest = load_manifest(site_manifests)
    site_folders = list(site_folders) if site_folders is not None else [DEFAULT_SITE_DIR]

    def site_file_for(species):
        """Path of a specimen source file: the first folder that has it, else where to save it (first folder)."""
        def f(g, lv, k):
            name = f"atlas_pair_{g}_{safe_site(lv)}_{k}.xlsx"
            for d in site_folders:
                p = os.path.join(d, safe(species), name)
                if os.path.exists(p):
                    return p
            return os.path.join(site_folders[0] if site_folders else DEFAULT_SITE_DIR, safe(species), name)
        return f
    downloads = {}   # follow-up files still needed: (kind, species, age group, level, file kind) -> info
    for r in rows:
        species, (a, b) = r["species"], r["comparison"].split("_vs_")
        bp = baseline_paths[species]
        if (species, a, b) not in panel_cache:
            panel_cache[(species, a, b)] = panel_coverage(bp[a], bp[b], PANEL_NAME.get(species, species), None)
            comp_cache[(species, a, b)] = composition_flags(bp[a], bp[b])
        pc, comp = panel_cache[(species, a, b)], comp_cache[(species, a, b)]
        r.update({k: v for k, v in pc.items() if not k.startswith("_")})
        r["panel_coverage_note"] = panel_note(pc, a, b, r["phenotype_flag"], r["country_adjusted_rate_diff_pp"])
        r.update(specimen_source_tvd_pp=comp["specimen_source"]["tvd_pp"],
                 specimen_source_differs=comp["specimen_source"]["differs"],
                 year_tvd_pp=comp["year"]["tvd_pp"], year_differs=comp["year"]["differs"])
        if r["status"] != "robust":
            continue
        kind = r["phenotype_flag"].lower()

        if not r["specimen_source_differs"]:
            r.update(specimen_source_check="not needed (specimen source mix matches between age groups)",
                     robust_after_specimen_source_check=True)
        else:
            site_file = site_file_for(species)
            sites = manifest.get(species) or {g: sorted({lv for d in site_folders
                                                         for lv in levels_from_files(os.path.join(d, safe(species)), g)})
                                              for g in (a, b)}
            cols, ok = stratified_rate_check(
                PANEL_NAME.get(species, species), kind, a, b, sites.get(a, []), sites.get(b, []), site_file)
            if ok == "unknown":
                missing = [x for x in needed_files("Specimen Source Group", kind, a, b, bp[a], bp[b], site_file) if not x[4]]
                for g, lv, k, p, _ in missing:
                    downloads.setdefault(("source", species, g, lv, k), {"path": p, "findings": set()})
                    downloads[("source", species, g, lv, k)]["findings"].add(f"{r['phenotype_flag']} {r['comparison']}")
                if missing:
                    cols["check"] = (f"flagged, files needed: download {len(missing)} specimen source file(s) listed in "
                                     f"files_to_download.txt, then run again")
            r.update(specimen_source_check=cols["check"].replace("adjustment", "specimen source adjustment"),
                     site_adjusted_rate_diff_pp=cols.get("adjusted_diff_pp", ""),
                     site_adjusted_ci_lo_pp=cols.get("ci_lo_pp", ""), site_adjusted_ci_hi_pp=cols.get("ci_hi_pp", ""),
                     n_sites_used=cols["n_used"], sites_used=cols.get("used", ""), per_site_diffs=cols.get("per_level", ""),
                     robust_after_specimen_source_check=ok)

        if not r["year_differs"]:
            r.update(year_check="not needed (year mix matches between age groups)", robust_after_year_check=True)
        else:
            yd = os.path.join(year_dir, safe(species))
            def years(g):
                pre = f"atlas_pair_{g}_"
                if not os.path.isdir(yd):
                    return []
                return sorted({f[len(pre):-len("_baseline.xlsx")] for f in os.listdir(yd)
                               if f.startswith(pre) and f.endswith("_baseline.xlsx")})
            year_file = lambda g, lv, k: os.path.join(yd, f"atlas_pair_{g}_{lv}_{k}.xlsx")
            cols, ok = stratified_rate_check(
                PANEL_NAME.get(species, species), kind, a, b, years(a), years(b), year_file)
            if ok == "unknown":
                missing = [x for x in needed_files("Year", kind, a, b, bp[a], bp[b], year_file) if not x[4]]
                for g, lv, k, p, _ in missing:
                    downloads.setdefault(("year", species, g, lv, k), {"path": p, "findings": set()})
                    downloads[("year", species, g, lv, k)]["findings"].add(f"{r['phenotype_flag']} {r['comparison']}")
                if missing:
                    cols["check"] = (f"flagged, files needed: download {len(missing)} year file(s) listed in "
                                     f"files_to_download.txt, then run again")
            r.update(year_check=cols["check"].replace("adjustment", "year adjustment"),
                     year_adjusted_rate_diff_pp=cols.get("adjusted_diff_pp", ""),
                     year_adjusted_ci_lo_pp=cols.get("ci_lo_pp", ""), year_adjusted_ci_hi_pp=cols.get("ci_hi_pp", ""),
                     n_years_used=cols["n_used"], years_used=cols.get("used", ""), per_year_diffs=cols.get("per_level", ""),
                     robust_after_year_check=ok)

        checks = (r["robust_after_specimen_source_check"], r["robust_after_year_check"])
        r["robust_final"] = (True if all(c is True for c in checks)
                             else False if any(c is False for c in checks) else "unknown")

    # ---- Write outputs
    os.makedirs(out_dir, exist_ok=True)
    for name, data in (("mdr_all_comparisons.csv", rows), ("mdr_trends.csv", robust)):
        with open(os.path.join(out_dir, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=ALL_FIELDS, extrasaction="ignore")
            w.writeheader()
            w.writerows(data)
    def tally(col, data):
        out = {}
        for r in data:
            out[str(r[col])] = out.get(str(r[col]), 0) + 1
        return out
    summary = {"n_species_with_data": len(baseline_paths), "species_without_mdr_dtr_exports": skipped,
               "n_computed": len(rows), "bh_family_m": len(fam), "status_counts": tally("status", rows),
               "n_robust": len(robust), "robust_by_flag": tally("phenotype_flag", robust),
               "robust_panel_comparable": sum(1 for r in robust if r["panel_coverage_verdict"] == "comparable"),
               "flagged_specimen_source": sum(1 for r in robust if r["specimen_source_differs"]),
               "flagged_year": sum(1 for r in robust if r["year_differs"]),
               "robust_after_specimen_source_check": tally("robust_after_specimen_source_check", robust),
               "robust_after_year_check": tally("robust_after_year_check", robust),
               "n_robust_final": tally("robust_final", robust).get("True", 0)}
    with open(os.path.join(out_dir, "mdr_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"Wrote mdr_all_comparisons.csv, mdr_trends.csv, mdr_summary.json to {out_dir}")

    # ---- Follow-up files still needed
    dl_path = os.path.join(out_dir, "files_to_download.txt")
    if not downloads:
        if os.path.exists(dl_path):
            os.remove(dl_path)
        return summary
    waiting = [r for r in robust if r["robust_final"] == "unknown"]
    lines = ["FILES TO DOWNLOAD FOR THE MDR/DTR SPECIMEN SOURCE / YEAR CHECKS", "",
             f"{len(waiting)} robust MDR/DTR finding(s) have a specimen source or year mix that differs between the age",
             "groups. They need the files below before they can be confirmed; until then their robust_final is 'unknown'.",
             "Download each file from ATLAS (MIC Distribution), save it with exactly the name and folder shown, then run",
             "mdr_full_pipeline.py again.", "",
             "How to download each file:",
             "  - Select the bacterium shown and ALL antibiotics, ALL years (except for year files).",
             "  - Set the age filter to the age group shown, and turn View By: Country OFF (worldwide).",
             "  - Specimen source files: select ONLY the specimen source shown. Year files: ONLY the year shown.",
             "  - Phenotype filter: none for _baseline files, MDR (CLSI) for _mdr files, DTR for _dtr files.",
             "  - Export, and save the .xlsx under the name below. In names, spaces are _ and '/' is -.", ""]
    for kind, title in (("source", "SPECIMEN SOURCE FILES"), ("year", "YEAR FILES")):
        items = sorted((k, v) for k, v in downloads.items() if k[0] == kind)
        if not items:
            continue
        lines.append(f"==== {title} ====")
        for species in sorted({k[1] for k, _ in items}):
            sp_items = [(k, v) for k, v in items if k[1] == species]
            finds = sorted(set().union(*[v["findings"] for _, v in sp_items]))
            lines.append(f"{species}  (needed for: {'; '.join(finds)})")
            for (_, _, g, lv, k), v in sp_items:
                lines.append(f"  age {g:6s} {lv:28s} {k:9s} -> {v['path']}")
            lines.append("")
    with open(dl_path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"\n>>> ACTION NEEDED: {len(waiting)} robust finding(s) need a specimen source or year check that this run could "
          f"not do, because {len(downloads)} file(s) are missing.\n>>> The list of files, where to save them and how to "
          f"download them is in {dl_path}\n>>> Download them, then run mdr_full_pipeline.py again.")
    return summary


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
    ap.add_argument("input", nargs="?", default=os.path.join(HERE, "mdr_input.txt"), help="input file (default: mdr_input.txt)")
    args = ap.parse_args()
    s = read_settings(args.input, need_bacteria=True)
    print(f"Input file: {s['path']}")
    run(species_list=sorted(s["pairs"]) if s["pairs"] is not None else None,
        data_dir=require_folder(s, "mdr_dtr_folder", "MDR/DTR"), site_folders=s["specimen_source_folders"] or None,
        year_dir=s["year_folder"] or os.path.join(os.path.dirname(s["path"]), "data", "year"),
        out_dir=s["results_folder"])


if __name__ == "__main__":
    main()
