"""
MRSA / CRAB marker pipeline: do single resistant lineages explain the S. aureus and A. baumannii age shifts?

Markers (the paper's definitions): MRSA = oxacillin MIC >= 4; CRAB = meropenem MIC >= 8, and separately
imipenem MIC >= 8. Assumes the three age-group exports (View By: Country, all drugs) are on disk.

Parts:
  1. Prevalence by age (S5 Table): country-adjusted marker prevalence in each age group and the difference,
     countries with >= 20 isolates of the marker drug in both age groups (>= 3 needed), 20,000-draw country
     bootstrap (binomial), BH across the 9 marker comparisons.
  2. Country-level correlation (Figures 6 and 7), 0-18 vs 65+: Spearman correlation across countries between the
     marker gap (older minus younger, percentage points) and each drug's log2 MIC shift, plus the mean shift over
     a drug group (S. aureus beta-lactams; A. baumannii non-carbapenem drugs). A country is used when it has >= 20
     isolates in both age groups for the marker drug and for the drug. Minimum countries: 10 (S. aureus),
     15 (A. baumannii, per-drug).
  3. Regression at zero marker gap (S9 Table): per drug and age comparison, weighted least squares of country
     shift on country marker gap (weight = harmonic-mean N); the intercept is the shift expected with no marker
     difference. 20,000 country bootstrap draws; >= 10 countries needed. Marker-defining drugs are left out.
  4. Direct removal (S10 Table): for the marker-defining drug only, drop isolates at or above the cutoff (MRSA or
     CRAB) and rerun the country-adjusted shift and bootstrap on what is left, over the same countries as the
     all-isolate shift.

Outputs, in results_folder:
  marker_prevalence.csv, marker_country_correlations.csv, marker_zero_gap_intercepts.csv, marker_removal.csv

Settings are in marker_input.txt. Run padi_full_pipeline.py first so the zero-gap table can mark which rows are
robust findings (from its padi_age_trends.csv).

Usage:
  python3 marker_pipeline.py                 (reads marker_input.txt)
  python3 marker_pipeline.py my_input.txt
"""
import argparse
import math
import os
import re

import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import spearmanr

# ============================================================================
# Settings (the paper's values)
# ============================================================================
MIN_COUNTRY_N = 20
MIN_QUALIFYING_COUNTRIES = 3
N_BOOTSTRAP = 2000
N_BOOTSTRAP_GATE = 20000
N_BOOT = 20000
MIN = MIN_COUNTRY_N
RNG_SEED = 0
REGRESSION_SEED = 20261114
AGE_GROUPS = ["0-18", "19-64", "65+"]
COMPARISONS = [("0-18", "19-64"), ("19-64", "65+"), ("0-18", "65+")]
COMPS = [("0-18", "19-64"), ("0-18", "65+"), ("19-64", "65+")]
TITLE_RE = re.compile(r"from (.+?)(?:\s*\(N|$)")

MARKERS = [   # part 1
    {"marker": "MRSA (oxacillin MIC >= 4)", "species": "Staphylococcus aureus", "drug": "Oxacillin", "threshold": 4.0},
    {"marker": "Carbapenem-resistant A. baumannii (meropenem MIC >= 8)", "species": "Acinetobacter baumannii",
     "drug": "Meropenem", "threshold": 8.0},
    {"marker": "Carbapenem-resistant A. baumannii (imipenem MIC >= 8)", "species": "Acinetobacter baumannii",
     "drug": "Imipenem", "threshold": 8.0},
]
CARBAPENEMS = ["Meropenem", "Imipenem", "Doripenem", "Ertapenem", "Tebipenem", "Meropenem-vaborbactam"]
SPECIES = {   # parts 3 and 4
    "Staphylococcus aureus": {"marker": "MRSA", "drug": "Oxacillin", "cut": 4, "exclude": ["Oxacillin"]},
    "Acinetobacter baumannii": {"marker": "CRAB (meropenem)", "drug": "Meropenem", "cut": 8, "exclude": CARBAPENEMS},
    "Acinetobacter baumannii (imipenem CRAB)": {"species": "Acinetobacter baumannii", "marker": "CRAB (imipenem)",
                                                "drug": "Imipenem", "cut": 8, "exclude": CARBAPENEMS},
}
SAUREUS_BETA_LACTAMS = ["Meropenem", "Doripenem", "Ertapenem", "Imipenem", "Ceftriaxone",
                        "Piperacillin-tazobactam", "Amoxicillin-clavulanate"]

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA_DIR = os.path.join(HERE, "data", "age_groups")
DEFAULT_OUT_DIR = os.path.join(HERE, "results", "markers")

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


def geomean_from_counts(counts):
    total = sum(counts.values())
    if total == 0:
        return None
    log2_sum = sum(math.log2(mic) * c for mic, c in counts.items())
    return 2 ** (log2_sum / total)


def weighted_shift(entries):
    shifts, weights = [], []
    for _, ea, eb in entries:
        ga, gb = geomean_from_counts(ea["counts"]), geomean_from_counts(eb["counts"])
        if not ga or not gb:
            continue
        shifts.append(math.log2(gb / ga))
        weights.append(2 * ea["N"] * eb["N"] / (ea["N"] + eb["N"]))
    if not shifts:
        return None
    return sum(w * s for w, s in zip(weights, shifts)) / sum(weights)


def bootstrap_p_value(boot, point):
    """Two-sided bootstrap p-value for H0: shift = 0, from the same draws that give the percentile CI.
    p = 2 * (1 + #draws on the far side of zero from the point estimate) / (B + 1), capped at 1.
    The +1 keeps p > 0 (a finite bootstrap can never prove p = 0). p < 0.05 matches 'CI excludes 0'
    up to that +1 correction."""
    boot = np.asarray(boot, dtype=float)
    if point is None or boot.size == 0:
        return None
    far = int((boot <= 0).sum()) if point > 0 else int((boot >= 0).sum()) if point < 0 else boot.size
    return min(1.0, 2.0 * (far + 1) / (boot.size + 1))


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


def _grid(ea, eb):
    """Common rung grid for one country's two age groups (sorted MIC values) + count vectors on it."""
    mics = sorted(set(ea["counts"]) | set(eb["counts"]))
    ca = np.array([ea["counts"].get(m, 0) for m in mics], dtype=float)
    cb = np.array([eb["counts"].get(m, 0) for m in mics], dtype=float)
    return np.array(mics, dtype=float), ca, cb


def _prob_superiority(ca, cb):
    """P(MIC_b > MIC_a) + 0.5 * P(tie), from count arrays on a common sorted grid (last axis).
    Rank-based, so it is unchanged by censoring at the panel edges (censored values keep their order;
    they only become ties). 0.5 = no difference; >0.5 = higher MICs in group b."""
    pa = ca / ca.sum(axis=-1, keepdims=True)
    pb = cb / cb.sum(axis=-1, keepdims=True)
    below_a = np.cumsum(pa, axis=-1) - pa          # P(MIC_a < rung)
    return (pb * (below_a + 0.5 * pa)).sum(axis=-1)


def weighted_prob_superiority(entries):
    """Country-weighted probability of superiority, same harmonic-mean-N weights as weighted_shift."""
    vals, weights = [], []
    for _, ea, eb in entries:
        _, ca, cb = _grid(ea, eb)
        if ca.sum() == 0 or cb.sum() == 0:
            continue
        vals.append(float(_prob_superiority(ca, cb)))
        weights.append(2 * ea["N"] * eb["N"] / (ea["N"] + eb["N"]))
    if not vals:
        return None
    return sum(w * v for w, v in zip(weights, vals)) / sum(weights)


def bootstrap_ci(entries, n_boot=N_BOOTSTRAP, seed=RNG_SEED, with_rank=False):
    """Two-level cluster bootstrap (unchanged design, vectorized 2026-11-08 so 20,000 draws are cheap):
    resample countries with replacement, then multinomial-resample each picked country's isolates
    (fixed N, observed MIC distribution); recompute the harmonic-mean-N-weighted shift per draw.
    Each pick of a country gets its own independent isolate resample, as in the old loop version.
    Returns the 95% percentile CI plus a two-sided bootstrap p-value (bootstrap_p_value).
    with_rank=True also returns the same statistics for the rank-based effect (probability of
    superiority, null 0.5), computed from the SAME resampled draws."""
    rng = np.random.default_rng(seed)
    n_c = len(entries)
    if n_c == 0:
        return None
    idx = rng.integers(0, n_c, size=(n_boot, n_c))
    mean_a = np.empty((n_boot, n_c))
    mean_b = np.empty((n_boot, n_c))
    ps = np.empty((n_boot, n_c)) if with_rank else None
    w_c = np.empty(n_c)
    for k, (_, ea, eb) in enumerate(entries):
        w_c[k] = 2 * ea["N"] * eb["N"] / (ea["N"] + eb["N"])
        mask = idx == k
        n_pick = int(mask.sum())
        if n_pick == 0:
            continue
        mics, ca, cb = _grid(ea, eb)
        logm = np.log2(mics)
        da = rng.multinomial(int(ea["N"]), ca / ca.sum(), size=n_pick)
        db = rng.multinomial(int(eb["N"]), cb / cb.sum(), size=n_pick)
        mean_a[mask] = (da @ logm) / da.sum(axis=1)
        mean_b[mask] = (db @ logm) / db.sum(axis=1)
        if with_rank:
            ps[mask] = _prob_superiority(da.astype(float), db.astype(float))
    w = w_c[idx]
    arr = (w * (mean_b - mean_a)).sum(axis=1) / w.sum(axis=1)
    lo, hi = np.percentile(arr, [2.5, 97.5])
    out = {"ci_lo": round(float(lo), 3), "ci_hi": round(float(hi), 3),
           "excludes_zero": bool(lo > 0 or hi < 0),
           "p_value": bootstrap_p_value(arr, weighted_shift(entries)), "n_boot": int(n_boot)}
    if with_rank:
        parr = (w * ps).sum(axis=1) / w.sum(axis=1)
        point = weighted_prob_superiority(entries)
        plo, phi = np.percentile(parr, [2.5, 97.5])
        out["rank"] = {"prob_superiority": round(point, 4), "ci_lo": round(float(plo), 4), "ci_hi": round(float(phi), 4),
                       "p_value": bootstrap_p_value(parr - 0.5, point - 0.5)}
    return out



# ---- rate difference and bootstrap (the same functions as in mdr_full_pipeline.py)
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



# ---- zero-gap regression and marker removal
def pct_ge(e, cut):
    return 100 * sum(n for m, n in e["counts"].items() if m >= cut) / e["N"]


def log2_shift(ea, eb):
    return math.log2(geomean_from_counts(eb["counts"]) / geomean_from_counts(ea["counts"]))


def wls(x, y, w):
    X = np.column_stack([np.ones_like(x), x])
    W = w / w.sum()
    beta = np.linalg.lstsq(X * np.sqrt(W)[:, None], y * np.sqrt(W), rcond=None)[0]
    return beta  # intercept, slope


def intercept_test(by, a, b, cfg, drug, rng):
    xs, ys, ws = [], [], []
    for c in sorted(set(by[a]) & set(by[b])):
        ma, mb = by[a][c].get(cfg["drug"]), by[b][c].get(cfg["drug"])
        da, db = by[a][c].get(drug), by[b][c].get(drug)
        if not (ma and mb and da and db) or min(ma["N"], mb["N"], da["N"], db["N"]) < MIN:
            continue
        xs.append(pct_ge(mb, cfg["cut"]) - pct_ge(ma, cfg["cut"]))
        ys.append(log2_shift(da, db))
        ws.append(2 * da["N"] * db["N"] / (da["N"] + db["N"]))
    n = len(xs)
    if n < 10:
        return None
    x, y, w = map(np.array, (xs, ys, ws))
    a0, b0 = wls(x, y, w)
    plain = float((w * y).sum() / w.sum())
    idx = rng.integers(0, n, size=(N_BOOT, n))
    boots = []
    for row in idx:
        xb = x[row]
        if np.ptp(xb) == 0:
            continue
        boots.append(wls(xb, y[row], w[row]))
    boots = np.array(boots)
    lo, hi = np.percentile(boots[:, 0], [2.5, 97.5])
    far = (boots[:, 0] <= 0).sum() if a0 > 0 else (boots[:, 0] >= 0).sum()
    return {"n_countries": n, "plain_weighted_shift": round(plain, 3), "intercept_at_zero_gap": round(float(a0), 3),
            "intercept_ci_lo": round(float(lo), 3), "intercept_ci_hi": round(float(hi), 3),
            "intercept_p": round(min(1.0, 2 * (far + 1) / (len(boots) + 1)), 5),
            "slope_per_10pp": round(float(b0) * 10, 3),
            "share_of_shift_left_at_zero_gap": round(float(a0) / plain, 2) if abs(plain) > 1e-9 else None,
            "mean_marker_gap_pp": round(float((w * x).sum() / w.sum()), 1)}


def stratum_test(by, a, b, cfg):
    def trunc(e):
        counts = {m: n for m, n in e["counts"].items() if m < cfg["cut"]}
        return {"N": sum(counts.values()), "counts": counts}
    full, strat = [], []
    for c in sorted(set(by[a]) & set(by[b])):
        ea, eb = by[a][c].get(cfg["drug"]), by[b][c].get(cfg["drug"])
        if not (ea and eb) or min(ea["N"], eb["N"]) < MIN:
            continue
        ta, tb = trunc(ea), trunc(eb)
        if min(ta["N"], tb["N"]) >= MIN:      # same countries for both rows, so "% of shift lost" compares like with like
            full.append((c, ea, eb))
            strat.append((c, ta, tb))
    out = {}
    for lab, ent in [("all_isolates", full), ("marker_free_only", strat)]:
        if len(ent) < MIN_QUALIFYING_COUNTRIES:
            out[lab] = None
            continue
        bt = bootstrap_ci(ent, n_boot=N_BOOT)
        lowest = min(m for _, ea, eb in ent for m in list(ea["counts"]) + list(eb["counts"]))
        at_low = sum(e["counts"].get(lowest, 0) for _, ea, eb in ent for e in (ea, eb))
        tot = sum(e["N"] for _, ea, eb in ent for e in (ea, eb))
        out[lab] = {"n_countries": len(ent), "n_isolates_a": sum(ea["N"] for _, ea, _ in ent),
                    "n_isolates_b": sum(eb["N"] for _, _, eb in ent),
                    "shift": round(weighted_shift(ent), 3), "ci_lo": bt["ci_lo"], "ci_hi": bt["ci_hi"],
                    "p": bt["p_value"], "pct_at_lowest_rung": round(100 * at_low / tot, 1), "lowest_rung": lowest}
    if out.get("all_isolates") and out.get("marker_free_only"):
        lost = 1 - out["marker_free_only"]["shift"] / out["all_isolates"]["shift"]
        out["marker_free_only"]["pct_shift_lost"] = round(100 * lost)
    return out



# ============================================================================
# Helpers
# ============================================================================
def safe(s):
    return s.replace(" ", "_").replace(",", "")


def load_species(data_dir, species, groups=AGE_GROUPS):
    d = os.path.join(data_dir, safe(species))
    return {g: parse_by_country_blocks(os.path.join(d, f"atlas_pair_{g}.xlsx")) for g in groups
            if os.path.exists(os.path.join(d, f"atlas_pair_{g}.xlsx"))}


def adj_rate(entries, i):
    """Country-adjusted rate of one age group: per-country rate weighted by the same harmonic-mean N as
    weighted_rate_diff, so adj_rate(b) - adj_rate(a) == the country-adjusted difference exactly.
    entries: (country, n_a, pos_a, n_b, pos_b); i = 2 for group a, 4 for group b."""
    w = [2 * e[1] * e[3] / (e[1] + e[3]) for e in entries]
    return sum(wi * e[i] / e[i - 1] for wi, e in zip(w, entries)) / sum(w)


# ============================================================================
# Part 1: prevalence by age
# ============================================================================
def marker_prevalence(data_dir):
    rows = []
    for mk in MARKERS:
        by_age = load_species(data_dir, mk["species"])
        rungs = sorted({m for bc in by_age.values() for dd in bc.values() for m in (dd.get(mk["drug"]) or {}).get("counts", {})})
        for a, b in COMPARISONS:
            if a not in by_age or b not in by_age:
                continue
            entries, pos, tot = [], {a: 0, b: 0}, {a: 0, b: 0}
            for c in sorted(set(by_age[a]) & set(by_age[b])):
                ea, eb = by_age[a][c].get(mk["drug"]), by_age[b][c].get(mk["drug"])
                if not ea or not eb or ea["N"] < MIN_COUNTRY_N or eb["N"] < MIN_COUNTRY_N:
                    continue
                pa = sum(n for m, n in ea["counts"].items() if m >= mk["threshold"])
                pb = sum(n for m, n in eb["counts"].items() if m >= mk["threshold"])
                entries.append((c, ea["N"], pa, eb["N"], pb))
                pos[a] += pa; pos[b] += pb; tot[a] += ea["N"]; tot[b] += eb["N"]
            row = {"marker": mk["marker"], "species": mk["species"], "comparison": f"{a}_vs_{b}",
                   "panel_rungs": "; ".join(f"{m:g}" for m in rungs), "n_countries": len(entries)}
            if len(entries) >= MIN_QUALIFYING_COUNTRIES:
                boot = bootstrap_rate_ci(entries, n_boot=N_BOOTSTRAP_GATE)
                row.update(pooled_rate_a_pct=round(100 * pos[a] / tot[a], 2), pooled_rate_b_pct=round(100 * pos[b] / tot[b], 2),
                           country_adjusted_rate_a_pct=round(100 * adj_rate(entries, 2), 2),
                           country_adjusted_rate_b_pct=round(100 * adj_rate(entries, 4), 2),
                           country_adjusted_rate_diff_pp=round(100 * weighted_rate_diff(entries), 2),
                           bootstrap_ci_lo_pp=boot["ci_lo_pp"], bootstrap_ci_hi_pp=boot["ci_hi_pp"],
                           bootstrap_p_value=boot["p_value"])
            rows.append(row)
    fam = [r for r in rows if r.get("bootstrap_p_value") is not None]
    for r, q in zip(fam, bh_qvalues([r["bootstrap_p_value"] for r in fam])):
        r["bh_q_value"] = q
    return pd.DataFrame(rows)


# ============================================================================
# Part 2: country-level correlation (0-18 vs 65+)
# ============================================================================
def country_correlations(data_dir):
    a, b = "0-18", "65+"
    out = []

    def table(by, marker_drug, cut, drugs, group_drugs, group_name):
        def ok(c, dr):
            ea, eb = by[a][c].get(dr), by[b][c].get(dr)
            return ea and eb and ea["N"] >= MIN_COUNTRY_N and eb["N"] >= MIN_COUNTRY_N
        def shift(c, dr):
            return math.log2(geomean_from_counts(by[b][c][dr]["counts"]) / geomean_from_counts(by[a][c][dr]["counts"]))
        def pct(e):
            return 100 * sum(n for m, n in e["counts"].items() if m >= cut) / e["N"]
        rows = []
        for c in sorted(set(by[a]) & set(by[b])):
            if not ok(c, marker_drug):
                continue
            row = {"gap": pct(by[b][c][marker_drug]) - pct(by[a][c][marker_drug])}
            g = [shift(c, d) for d in group_drugs if ok(c, d)]
            row[group_name] = np.mean(g) if g else np.nan
            for d in drugs:
                row[d] = shift(c, d) if ok(c, d) else np.nan
            rows.append(row)
        return pd.DataFrame(rows)

    def corr(df, col, species, marker, min_n, kind):
        x = df.dropna(subset=[col])
        if len(x) < min_n:
            return
        rho, p = spearmanr(x["gap"], x[col])
        out.append({"species": species, "marker": marker, "comparison": f"{a}_vs_{b}", "drug_or_group": col,
                    "kind": kind, "rho": round(float(rho), 3), "p_value": float(p), "n_countries": len(x)})

    # S. aureus: MRSA gap vs beta-lactam mean and each drug
    by = load_species(data_dir, "Staphylococcus aureus", (a, b))
    drugs = sorted({d for c in by[a].values() for d in c} - {"Oxacillin"})
    t = table(by, "Oxacillin", 4, drugs, SAUREUS_BETA_LACTAMS, "beta-lactam mean shift")
    corr(t, "beta-lactam mean shift", "Staphylococcus aureus", "MRSA", 10, "group mean")
    for d in drugs:
        corr(t, d, "Staphylococcus aureus", "MRSA", 10, "drug")

    # A. baumannii: CRAB gap (meropenem, imipenem) vs non-carbapenem mean and each drug
    by = load_species(data_dir, "Acinetobacter baumannii", (a, b))
    drugs = sorted({d for c in by[a].values() for d in c})
    group = [d for d in drugs if d not in CARBAPENEMS + ["Colistin"]]
    for crab_drug in ("Meropenem", "Imipenem"):
        t = table(by, crab_drug, 8, drugs, group, "non-carbapenem mean shift")
        marker = f"CRAB ({crab_drug.lower()})"
        corr(t, "non-carbapenem mean shift", "Acinetobacter baumannii", marker, 10, "group mean")
        for d in drugs:
            corr(t, d, "Acinetobacter baumannii", marker, 15, "drug")
    return pd.DataFrame(out)


# ============================================================================
# Parts 3 and 4: regression at zero marker gap, direct removal
# ============================================================================
def zero_gap_and_removal(data_dir):
    rng = np.random.default_rng(REGRESSION_SEED)
    ints, strata = [], []
    for key, cfg in SPECIES.items():
        sp = cfg.get("species", key)
        by = load_species(data_dir, sp)
        drugs = sorted({x for bc in by.values() for dd in bc.values() for x in dd} - set(cfg["exclude"]))
        for a, b in COMPS:
            for drug in drugs:
                r = intercept_test(by, a, b, cfg, drug, rng)
                if r:
                    ints.append({"species": sp, "marker": cfg["marker"], "comparison": f"{a}_vs_{b}",
                                 "antibacterial": drug, **r})
            s = stratum_test(by, a, b, cfg)
            for lab, v in s.items():
                if v:
                    strata.append({"species": sp, "marker": cfg["marker"], "marker_drug": cfg["drug"],
                                   "comparison": f"{a}_vs_{b}", "isolates": lab, **v})
    return pd.DataFrame(ints), pd.DataFrame(strata)


def robust_keys(path):
    """(species, drug, comparison) of final robust findings, from padi_full_pipeline.py's padi_age_trends.csv
    (column robust_final) or the older robust_signals_PADI.csv (robust_after_specimen_source_adjustment)."""
    if not path or not os.path.exists(path):
        return None
    r = pd.read_csv(path)
    col = "robust_final" if "robust_final" in r.columns else "robust_after_specimen_source_adjustment"
    r = r[r[col].astype(str) == "True"]
    return set(zip(r.species, r.antibacterial, r.comparison))


def run(data_dir=DEFAULT_DATA_DIR, out_dir=DEFAULT_OUT_DIR, robust_table=None):
    missing = [os.path.join(data_dir, safe(sp), f"atlas_pair_{g}.xlsx")
               for sp in ("Staphylococcus aureus", "Acinetobacter baumannii") for g in AGE_GROUPS
               if not os.path.exists(os.path.join(data_dir, safe(sp), f"atlas_pair_{g}.xlsx"))]
    if missing:
        raise SystemExit("This step needs the age group exports for S. aureus and A. baumannii. Missing:\n  "
                         + "\n  ".join(missing) + "\nSee marker_input.txt for how to download and name them.")
    os.makedirs(out_dir, exist_ok=True)
    print("Part 1: marker prevalence by age", flush=True)
    prev = marker_prevalence(data_dir)
    prev.to_csv(os.path.join(out_dir, "marker_prevalence.csv"), index=False)
    print(prev[["marker", "comparison", "country_adjusted_rate_a_pct", "country_adjusted_rate_b_pct",
                "country_adjusted_rate_diff_pp", "bootstrap_ci_lo_pp", "bootstrap_ci_hi_pp"]].to_string(index=False))

    print("\nPart 2: country-level correlations (0-18 vs 65+)", flush=True)
    cor = country_correlations(data_dir)
    cor.to_csv(os.path.join(out_dir, "marker_country_correlations.csv"), index=False)
    print(cor[cor.kind == "group mean"].to_string(index=False))

    print("\nParts 3-4: regression at zero marker gap, direct removal", flush=True)
    ints, strata = zero_gap_and_removal(data_dir)
    keys = robust_keys(robust_table)
    if keys is not None:
        ints["robust_finding"] = [(r.species, r.antibacterial, r.comparison) in keys for r in ints.itertuples()]
    ints.to_csv(os.path.join(out_dir, "marker_zero_gap_intercepts.csv"), index=False)
    strata.to_csv(os.path.join(out_dir, "marker_removal.csv"), index=False)
    print(strata[["marker", "comparison", "isolates", "shift", "ci_lo", "ci_hi"]
                 + (["pct_shift_lost"] if "pct_shift_lost" in strata else [])].to_string(index=False))
    if keys is not None:
        print(f"{int(ints.robust_finding.sum())} robust findings assessed by the zero-gap regression")
    print(f"\nWrote marker_prevalence.csv, marker_country_correlations.csv, marker_zero_gap_intercepts.csv, "
          f"marker_removal.csv to {out_dir}")


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
    ap.add_argument("input", nargs="?", default=os.path.join(HERE, "marker_input.txt"), help="input file (default: marker_input.txt)")
    args = ap.parse_args()
    s = read_settings(args.input)
    print(f"Input file: {s['path']}")
    robust = s["padi_age_trends"]
    if not robust or not os.path.exists(robust):
        print(f"NOTE: padi_age_trends not found ({robust}); run padi_full_pipeline.py first to mark robust findings.")
    run(require_folder(s, "age_group_folder", "age group"), s["results_folder"], robust)


if __name__ == "__main__":
    main()
