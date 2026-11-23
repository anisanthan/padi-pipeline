"""
PADI full pipeline: every step of the paper's MIC age analysis in one run.

Give it species (and optionally drugs). Assuming the ATLAS exports are already on disk, it outputs
the age trends: every pediatric-vs-adult comparison tested and which ones are robust.

Steps, in order (same math and cutoffs as the paper):
  1. Pediatric sample size check: keep only species-drug pairs marked "viable" in the sample size table
     (>= 20 pediatric isolates in >= 3 countries). Skip with --no-sample-size-filter.
  2. Country-adjusted log2 MIC shift for every pair and age comparison (0-18 vs 19-64, 0-18 vs 65+,
     19-64 vs 65+), using countries with >= 20 isolates in both age groups; >= 3 such countries needed.
  3. Country bootstrap (20,000 draws, fixed seed) -> 95% CI and two-sided p-value, for every comparison.
  4. Benjamini-Hochberg correction across every comparison tested in this run (one family).
     Robust = |shift| >= 0.5 log2 AND q < 0.05.
  5. Mixed effects model (reported, not gating) for every comparison that clears the 0.5 floor.
  6. For candidates (|shift| >= 0.5 and uncorrected p < 0.05, i.e. robust plus those that only fail BH):
     specimen source and year composition check, % resistant in each age group (country-adjusted), the
     EUCAST expected-resistant flag, and the median MIC of each age group with its susceptible / intermediate /
     resistant class (ATLAS cell colours). A factor is flagged when its mix differs between the two age groups by
     TVD > 15 pp with chi-square p < 0.01.
  7. Specimen source check for flagged candidates: shift within each specimen source, weighted and bootstrapped.
     Year check for flagged candidates: shift within each year (same method).
     For both checks, the finding survives only if the adjusted CI excludes zero; with < 3 usable specimen
     sources or years it is "not checked". robust_final = robust AND survives (or did not need) both checks.
     Calibration of the 15 pp threshold: across every candidate re-tested within specimen sources, the largest
     |pooled shift - source-adjusted shift| per 15 pp of specimen source difference.
  8. Specimen source mix (% of isolates from each source) per age group, for species with a flagged candidate.
  9. Optional permutation negative control (--permutations N): the same null datasets as the paper's
     control, run over the whole BH family; reports each comparison's null rate of |shift| >= 0.5 and p < 0.05.

The BH family is every comparison tested in THIS run, so q-values (and the robust set) match the paper
only when the run covers the same pairs (the default: every species folder in the data directory, all drugs).

Inputs: set in padi_input.txt (age_group_folder, sample_size_table from pediatric_sample_size_check.py,
specimen_source_folders, year_folder, results_folder), with the bacteria and drugs to analyze.
Outputs, in results_folder:
  padi_all_comparisons.csv   every comparison with its stage: insufficient_countries / below_effect_floor /
                             not_significant / fails_bh_fdr / robust
  padi_age_trends.csv        robust comparisons with direction, shift, CI, q, % resistant, median MIC class,
                             EUCAST flag, specimen source and year results, robust_final
  padi_specimen_source_mix.csv   specimen source percentages per age group (species with a flagged candidate)
  padi_summary.json          counts at each step, the 15 pp calibration and the median MIC class changes
  files_to_download.txt      only if a flagged finding still needs specimen source or year files

Usage:
  1. Run pediatric_sample_size_check.py first.
  2. In padi_input.txt (next to this script) set where your ATLAS files are and list the bacteria and drugs.
  3. Run:  python3 padi_full_pipeline.py                 (reads padi_input.txt)
     or:   python3 padi_full_pipeline.py my_input.txt
"""
import argparse
import csv
import itertools
import json
import math
import os
import re
import time
from multiprocessing import Pool

import numpy as np
import openpyxl
from scipy.stats import chi2_contingency

try:
    import statsmodels.formula.api as smf
    import pandas as pd
    HAVE_MIXEDLM = True
except ImportError:
    HAVE_MIXEDLM = False

# ============================================================================
# Settings (the paper's values)
# ============================================================================
MIN_GLOBAL_N = 20                 # isolates per age group per specimen source (specimen source check)
MIN_COUNTRY_N = 20                # isolates per age group per country
MIN_QUALIFYING_COUNTRIES = 3      # countries (or specimen sources) needed
EFFECT_FLOOR = 0.5                # |log2 shift| floor
N_BOOTSTRAP = 2000                # specimen source check bootstrap
N_BOOTSTRAP_GATE = 20000          # main country bootstrap
FDR_Q = 0.05                      # Benjamini-Hochberg level
RNG_SEED = 0
PERM_SEED = 20261114              # permutation control seed (same as the paper's control)
AGE_GROUPS = ["0-18", "19-64", "65+"]
COMPOSITION_SECTIONS = ["Specimen Source Group", "Year"]

COLOR_S = "FFA2F5B2"              # ATLAS cell fill colours for susceptible / intermediate / resistant
COLOR_I = "FFFFDA51"
COLOR_R = "FFFF8380"
TITLE_RE = re.compile(r"from (.+?)(?:\s*\(N|$)")

# ATLAS export drug names -> names in the sample size table
EXPORT_TO_TABLE_NAME = {
    "Trimethoprim-sulfamethoxazole": "Trimethoprim sulfa",
    "Amoxicillin-clavulanate": "Amoxy/clav",
    "Piperacillin-tazobactam": "Pip/taz",
    "Ampicillin-sulbactam": "Ampicillin sulbactam",
    "Cefoperazone-sulbactam": "Cefoperazone sulbactam",
    "Ceftazidime-avibactam": "Ceftazidime avibactam",
    "Ceftolozane-tazobactam": "Ceftolozane tazobactam",
    "Aztreonam-avibactam": "Aztreonam avibactam",
    "Meropenem-vaborbactam": "Meropenem vaborbactam",
    "Quinupristin-dalfopristin": "Quinupristin dalfopristin",
}

# EUCAST Expected Resistant Phenotypes v1.2 (January 2023), Tables 1-5, mapped to the drugs in the ATLAS exports.
# Used only to flag rows (expected_resistant_eucast); nothing is gated on it.
# https://www.eucast.org/fileadmin/eucast/pdf/expert_rules/Expected_Resistant_Phenotypes_v1.2_20230113.pdf
# Mapping choices (EUCAST columns -> the 48 drug names in our ATLAS exports):
#   - Only drugs present in our exports matter; EUCAST agents we never see (cefazolin, ticarcillin, trimethoprim
#     alone, fosfomycin, nitrofurantoin, chloramphenicol, fusidic acid, rifampicin, ...) are omitted.
#   - "Cefotaxime/Ceftriaxone" -> Ceftriaxone only (not extended to other 3rd-gen cephalosporins, e.g. cefixime).
#   - "Tetracyclines" -> Minocycline, EXCEPT Acinetobacter and S. maltophilia, where the footnotes say only
#     tetracycline/doxycycline are expected resistant (neither is in our exports).
#   - Class words in table headers -> our drugs: benzylpenicillin = Penicillin; glycopeptides = Vancomycin,
#     Teicoplanin; macrolides = Azithromycin, Erythromycin, Clarithromycin; lincosamides = Clindamycin;
#     oxazolidinones = Linezolid; aminoglycosides = Amikacin, Gentamicin, Tobramycin; first/second-generation
#     cephalosporins = Cefoxitin.
#   - Gram-positives: a beta-lactamase-inhibitor partner does not add gram-positive activity, so Aztreonam ->
#     also Aztreonam-avibactam; Ceftazidime -> also Ceftazidime-avibactam; "Cephalosporins" (enterococci) -> every
#     cephalosporin in our exports incl. inhibitor combinations and cefiderocol.
#   - Streptococci and enterococci aminoglycosides = LOW-LEVEL resistance (EUCAST footnote 1); still flagged.
#   - Ampicillin-sulbactam for A. baumannii is a "Note" (may test susceptible), NOT R -> not flagged.
#   - Trimethoprim-sulfamethoxazole is never flagged: EUCAST lists trimethoprim alone (P. aeruginosa,
#     Acinetobacter, S. maltophilia, M. catarrhalis) or sulfonamides alone (enterococci), not the combination.
PEN = {"Penicillin"}
GLYCO = {"Vancomycin", "Teicoplanin"}
MACRO = {"Azithromycin", "Erythromycin", "Clarithromycin"}
LINC = {"Clindamycin"}
OXAZ = {"Linezolid"}
AMINO = {"Amikacin", "Gentamicin", "Tobramycin"}
CEPH_1_2 = {"Cefoxitin"}
ALL_CEPHS = {"Cefepime", "Cefiderocol", "Cefixime", "Cefoperazone-sulbactam", "Cefoxitin", "Cefpodoxime",
             "Ceftaroline", "Ceftaroline-avibactam", "Ceftazidime", "Ceftazidime-avibactam", "Ceftibuten",
             "Ceftibuten-avibactam", "Ceftolozane-tazobactam", "Ceftriaxone"}

ENTEROBACTERALES = PEN | GLYCO | MACRO | LINC | OXAZ                      # Table 1 header
NONFERMENTER = PEN | CEPH_1_2 | GLYCO | MACRO | LINC | OXAZ               # Table 2 header
OTHER_GN = GLYCO | LINC | OXAZ                                            # Table 3 header
GRAM_POS = {"Aztreonam", "Aztreonam-avibactam", "Colistin"}               # Table 4 header
ANAEROBE = {"Aztreonam", "Aztreonam-avibactam", "Colistin"} | AMINO      # Table 5 header
AMPC = {"Ampicillin", "Amoxicillin-clavulanate", "Ampicillin-sulbactam", "Cefoxitin"}   # rules 1.2, 1.3, 1.6
CAZ = {"Ceftazidime", "Ceftazidime-avibactam"}

EXPECTED_RESISTANT = {
    # Table 1
    "Escherichia coli": ENTEROBACTERALES,
    "Citrobacter koseri": ENTEROBACTERALES | {"Ampicillin"},                                   # 1.1
    "Citrobacter freundii": ENTEROBACTERALES | AMPC,                                           # 1.2
    "Enterobacter cloacae": ENTEROBACTERALES | AMPC,                                           # 1.3 complex
    "Enterobacter hormaechei": ENTEROBACTERALES | AMPC,                                        # 1.3 complex
    "Klebsiella aerogenes": ENTEROBACTERALES | AMPC,                                           # 1.6
    "Klebsiella pneumoniae": ENTEROBACTERALES | {"Ampicillin"},                                # 1.7 complex
    "Klebsiella oxytoca": ENTEROBACTERALES | {"Ampicillin"},                                   # 1.8
    "Morganella morganii": ENTEROBACTERALES | {"Ampicillin", "Amoxicillin-clavulanate",
                                               "Ampicillin-sulbactam", "Minocycline", "Colistin"},   # 1.10
    "Proteus mirabilis": ENTEROBACTERALES | {"Minocycline", "Colistin"},                       # 1.12
    "Serratia marcescens": ENTEROBACTERALES | AMPC | {"Colistin"},                             # 1.18
    # Table 2
    "Acinetobacter baumannii": NONFERMENTER | {"Ampicillin", "Amoxicillin-clavulanate", "Ceftriaxone",
                                               "Aztreonam", "Ertapenem"},                       # 2.1
    "Acinetobacter pitii": NONFERMENTER | {"Ampicillin", "Amoxicillin-clavulanate", "Ceftriaxone",
                                           "Aztreonam", "Ertapenem"},                           # 2.1 (A. pittii)
    "Acinetobacter lwoffii": NONFERMENTER,                                                     # header only
    "Acinetobacter ursingii": NONFERMENTER,                                                    # header only
    "Pseudomonas aeruginosa": NONFERMENTER | {"Ampicillin", "Amoxicillin-clavulanate", "Ampicillin-sulbactam",
                                              "Ceftriaxone", "Ertapenem", "Minocycline", "Tigecycline"},  # 2.7
    "Stenotrophomonas maltophilia": NONFERMENTER | AMINO | {
        "Ampicillin", "Amoxicillin-clavulanate", "Ampicillin-sulbactam", "Piperacillin-tazobactam",
        "Ceftriaxone", "Aztreonam", "Ertapenem", "Imipenem", "Meropenem"},                     # 2.8
    # Table 3
    "Haemophilus influenzae": OTHER_GN,                                                        # 3.1
    "Moraxella catarrhalis": OTHER_GN,                                                         # 3.2
    # Table 4
    "Staphylococcus aureus": GRAM_POS | CAZ,                                                   # 4.5
    "Staphylococcus epidermidis": GRAM_POS | CAZ,                                              # 4.5
    "Staphylococcus haemolyticus": GRAM_POS | CAZ,                                             # 4.5
    "Streptococcus pneumoniae": GRAM_POS | CAZ | AMINO,                                        # 4.6 (LLR)
    "Streptococcus agalactiae": GRAM_POS | CAZ | AMINO,                                        # 4.6 (LLR)
    "Streptococcus pyogenes": GRAM_POS | CAZ | AMINO,                                          # 4.6 (LLR)
    "Enterococcus faecalis": GRAM_POS | ALL_CEPHS | AMINO | MACRO | LINC,                      # 4.7 (AMINO = LLR)
    "Enterococcus faecium": GRAM_POS | ALL_CEPHS | AMINO | MACRO,                              # 4.9 (AMINO = LLR)
    # Table 5
    "Bacteroides fragilis": ANAEROBE,
    "Clostridioides difficile": ANAEROBE,
}

LOW_LEVEL_ONLY = {(sp, d) for sp in ("Streptococcus pneumoniae", "Streptococcus agalactiae", "Streptococcus pyogenes",
                                      "Enterococcus faecalis", "Enterococcus faecium") for d in AMINO}


def is_expected_resistant(species, drug):
    return drug in EXPECTED_RESISTANT.get(species, set())


HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA_DIR = os.path.join(HERE, "data", "age_groups")
DEFAULT_SAMPLE_SIZE_TABLE = os.path.join(HERE, "results", "sample_size", "pediatric_sample_size.csv")
DEFAULT_SITE_DIRS = [   # searched in this order
    (os.path.join(HERE, "data", "specimen_source"), ""),
    (os.path.join(HERE, "data", "specimen_source_extra"), ""),
]
DEFAULT_SITE_MANIFESTS = []
DEFAULT_YEAR_DIR = os.path.join(HERE, "data", "year")
DEFAULT_OUT_DIR = os.path.join(HERE, "results", "padi")


# ---- parsing
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


def parse_by_country_breakpoints(path, sheet="Absolute"):
    """Same blocks, but reads cell fill color instead of value -- gives
    {country: {drug: {"S": n, "I": n, "R": n, "no_bp": n}}}."""
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(min_row=1, max_row=ws.max_row, values_only=False))
    data, current, header, in_table = {}, None, None, False
    for row in rows:
        first = row[0].value
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
        n_val = row[1].value
        if not isinstance(first, str) or not isinstance(n_val, (int, float)):
            in_table = False
            continue
        counts = {"S": 0, "I": 0, "R": 0, "no_bp": 0}
        for cell in row[2:]:
            v = cell.value
            if v is None or v == 0:
                continue
            fill = cell.fill.fgColor.rgb if cell.fill and cell.fill.patternType else None
            fill = fill.upper() if isinstance(fill, str) else fill
            key = {COLOR_S: "S", COLOR_I: "I", COLOR_R: "R"}.get(fill, "no_bp")
            counts[key] += int(v)
        data[current][first] = counts
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


def parse_flat_absolute(path, sheet="Absolute"):
    """Single-group worldwide export (no country/title blocks) --
    {drug: {"N": int, "counts": {mic: count}}}."""
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(min_row=1, max_row=ws.max_row, values_only=True))
    header, data = None, {}
    for row in rows:
        if row[0] == "Antibacterial":
            header = row
            continue
        if header is None or not isinstance(row[0], str) or not isinstance(row[1], (int, float)):
            continue
        counts = {}
        for mic_val, c in zip(header[2:], row[2:]):
            if mic_val is None or c is None:
                continue
            counts[float(mic_val)] = int(c)
        data[row[0]] = {"N": int(row[1]), "counts": counts}
    return data


# ---- math
def geomean_from_counts(counts):
    total = sum(counts.values())
    if total == 0:
        return None
    log2_sum = sum(math.log2(mic) * c for mic, c in counts.items())
    return 2 ** (log2_sum / total)


def multinomial_resample(counts, n, rng):
    mics = list(counts.keys())
    probs = np.array([counts[m] for m in mics], dtype=float)
    probs = probs / probs.sum()
    draws = rng.multinomial(n, probs)
    return {m: int(d) for m, d in zip(mics, draws)}


def qualifying_countries(by_a, by_b, drug):
    out = []
    for c in sorted(set(by_a) & set(by_b)):
        ea, eb = by_a[c].get(drug), by_b[c].get(drug)
        if not ea or not eb or ea["N"] < MIN_COUNTRY_N or eb["N"] < MIN_COUNTRY_N:
            continue
        out.append((c, ea, eb))
    return out


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


def fit_mixedlm(entries, ref_label, alt_label):
    if not HAVE_MIXEDLM:
        return {"converged": False, "reason": "statsmodels/pandas not installed"}
    rows = []
    for c, e_ref, e_alt in entries:
        for mic, cnt in e_ref["counts"].items():
            rows.append((c, ref_label, math.log2(mic), cnt))
        for mic, cnt in e_alt["counts"].items():
            rows.append((c, alt_label, math.log2(mic), cnt))
    df = pd.DataFrame(rows, columns=["country", "age_group", "log2_mic", "weight"])
    df = df.loc[df.index.repeat(df["weight"])][["country", "age_group", "log2_mic"]].reset_index(drop=True)
    if df.empty or df["country"].nunique() < 3:
        return {"converged": False, "reason": "too few countries with data"}
    df["age_group"] = pd.Categorical(df["age_group"], categories=[ref_label, alt_label])
    try:
        res = smf.mixedlm("log2_mic ~ age_group", df, groups=df["country"], re_formula="~age_group").fit(reml=True)
    except Exception:
        try:
            res = smf.mixedlm("log2_mic ~ age_group", df, groups=df["country"]).fit(reml=True)
        except Exception as e2:
            return {"converged": False, "reason": str(e2)}
    key = f"age_group[T.{alt_label}]"
    if key not in res.params.index:
        return {"converged": False, "reason": "coefficient not found in fit"}
    p = res.pvalues[key]
    return {"converged": True, "coef_log2_shift": round(float(res.params[key]), 3),
            "p_value": float(p), "excludes_zero": bool(p < 0.05)}


def make_null_entries(entries, rng):
    null_entries = []
    for country, ea, eb in entries:
        pooled = dict(ea["counts"])
        for mic, c in eb["counts"].items():
            pooled[mic] = pooled.get(mic, 0) + c
        mics = list(pooled.keys())
        probs = np.array([pooled[m] for m in mics], dtype=float)
        probs /= probs.sum()
        n_a, n_b = ea["N"], eb["N"]
        draw_a = dict(zip(mics, rng.multinomial(n_a, probs)))
        draw_b = dict(zip(mics, rng.multinomial(n_b, probs)))
        null_entries.append((country, {"N": n_a, "counts": draw_a}, {"N": n_b, "counts": draw_b}))
    return null_entries


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


# ---- site
def weighted_site_shift(entries):
    """Same math as weighted_shift, with specimen sources (or years) in place of countries:
    per-level log2 shift, weighted by the harmonic mean of the two group's
    isolate counts, averaged across every qualifying level."""
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


def bootstrap_site_ci(entries, n_boot=N_BOOTSTRAP, seed=RNG_SEED):
    """Same resampling scheme as bootstrap_ci, with specimen sources (or years) in place of countries:
    resample which levels are included (with replacement), multinomial-resample each level's isolate counts,
    recompute the weighted shift, repeat, take the middle 95% as the CI."""
    rng = np.random.default_rng(seed)
    n = len(entries)
    if n == 0:
        return None
    boot_shifts = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        resampled = []
        for i in idx:
            site, ea, eb = entries[i]
            ra = {"N": ea["N"], "counts": multinomial_resample(ea["counts"], ea["N"], rng)}
            rb = {"N": eb["N"], "counts": multinomial_resample(eb["counts"], eb["N"], rng)}
            resampled.append((site, ra, rb))
        s = weighted_site_shift(resampled)
        if s is not None:
            boot_shifts.append(s)
    if not boot_shifts:
        return None
    arr = np.array(boot_shifts)
    lo, hi = np.percentile(arr, [2.5, 97.5])
    return {"ci_lo": round(float(lo), 3), "ci_hi": round(float(hi), 3), "excludes_zero": bool(lo > 0 or hi < 0)}



# ============================================================================
# Helpers
# ============================================================================
def safe(s):
    return s.replace(" ", "_").replace(",", "")


def safe_site(s):
    return s.replace(" ", "_").replace("/", "-")


def norm_name(d):
    return "".join(ch for ch in d.lower() if ch.isalnum())


def load_viable_pairs(path):
    """{(species, normalized drug name)} for pairs that passed the pediatric sample size check."""
    rows = list(csv.DictReader(open(path)))
    return {(r["species"], norm_name(r["antibacterial"])) for r in rows if r["status"] == "viable"}


def is_viable(viable, species, drug):
    return (species, norm_name(EXPORT_TO_TABLE_NAME.get(drug, drug))) in viable


def load_site_manifest(paths):
    """{species: {age group: [specimen sources downloaded]}}, merged across manifest files."""
    m = {}
    for p in paths:
        if not os.path.exists(p):
            continue
        for sp, groups in json.load(open(p)).items():
            for g, sites in groups.items():
                have = m.setdefault(sp, {}).setdefault(g, [])
                have.extend(x for x in sites if x not in have)
    return m


def site_search_order(folders):
    """Every specimen source folder is searched for <age group>_<source>.xlsx first, then for
    <age group>_<source>_baseline.xlsx (the unfiltered export from the MDR/DTR specimen source downloads, which works
    the same way)."""
    return [(d, suf) for suf in ("", "_baseline") for d in folders]


def discover_sources(species, a, b, folders):
    """Specimen sources with a file for BOTH age groups, read from the file names
    <folder>/<Species>/atlas_pair_<age group>_<source>.xlsx (or ..._baseline.xlsx)."""
    def sources(g):
        out = set()
        for d in folders:
            sd = os.path.join(d, safe(species))
            if not os.path.isdir(sd):
                continue
            pre = f"atlas_pair_{g}_"
            for fn in os.listdir(sd):
                if fn.startswith(pre) and fn.endswith(".xlsx"):
                    name = fn[len(pre):-5]
                    if name.endswith(("_mdr", "_dtr")):
                        continue
                    name = name[:-len("_baseline")] if name.endswith("_baseline") else name
                    out.add(name.replace("_-_", " / ").replace("_", " "))   # file name -> source name
        return out
    return sorted(sources(a) & sources(b))


def gather_site_entries(species, drug, a, b, sites, site_dirs):
    """[(source, entry_a, entry_b)] for every specimen source with >= MIN_GLOBAL_N isolates of this drug in both
    age groups. Each source is looked up in site_dirs in order; the first folder with usable data is used."""
    entries = []
    for site in sites:
        ea = eb = None
        for d, suf in site_dirs:
            path_a = os.path.join(d, safe(species), f"atlas_pair_{a}_{safe_site(site)}{suf}.xlsx")
            path_b = os.path.join(d, safe(species), f"atlas_pair_{b}_{safe_site(site)}{suf}.xlsx")
            if not (os.path.exists(path_a) and os.path.exists(path_b)):
                continue
            ea = parse_flat_absolute(path_a).get(drug)
            eb = parse_flat_absolute(path_b).get(drug)
            if ea and eb and ea["N"] >= MIN_GLOBAL_N and eb["N"] >= MIN_GLOBAL_N:
                break
            ea = eb = None
        if ea and eb:
            entries.append((site, ea, eb))
    return entries


def country_adjusted_pct_resistant(blocks_a, blocks_b, bp_a, bp_b, drug):
    """% resistant in each age group over the comparison's qualifying countries, each country weighted by the
    harmonic mean of its two breakpoint-tested isolate counts. Returns (pct_a, pct_b, n_a, n_b) or None
    (no clinical breakpoint for this pair)."""
    num_a = num_b = wsum = 0.0
    n_a = n_b = 0
    for c, _, _ in qualifying_countries(blocks_a, blocks_b, drug):
        ea, eb = bp_a.get(c, {}).get(drug), bp_b.get(c, {}).get(drug)
        if not ea or not eb:
            continue
        da, db = ea["S"] + ea["I"] + ea["R"], eb["S"] + eb["I"] + eb["R"]
        if not da or not db:
            continue
        w = 2 * da * db / (da + db)
        num_a += w * ea["R"] / da
        num_b += w * eb["R"] / db
        wsum += w
        n_a += da
        n_b += db
    if not wsum:
        return None
    return round(100 * num_a / wsum, 1), round(100 * num_b / wsum, 1), n_a, n_b


def composition_flags(comp_a, comp_b):
    """Specimen source and year mix compared between two age groups (pooled over every country in the files).
    differs = TVD > 15 pp and chi-square p < 0.01."""
    out = {}
    for section in COMPOSITION_SECTIONS:
        cmp = compare_composition(pooled_category_counts(comp_a, section), pooled_category_counts(comp_b, section))
        if cmp is None:
            continue
        cmp["differs"] = bool(cmp["chi2_p"] is not None and cmp["tvd_pp"] > 15 and cmp["chi2_p"] < 0.01)
        out[section] = cmp
    return out


def discover_years(species, a, b, year_dir):
    """Year labels with a worldwide export for BOTH age groups:
    <year_dir>/<Species>/atlas_pair_<age group>_<year>.xlsx  (<year> = one year, e.g. 2019, or a range, e.g. 2004-2010)."""
    d = os.path.join(year_dir, safe(species))
    if not os.path.isdir(d):
        return []
    def labels(g):
        pre, out = f"atlas_pair_{g}_", set()
        for fn in os.listdir(d):
            if not (fn.startswith(pre) and fn.endswith(".xlsx")):
                continue
            name = fn[len(pre):-5]
            if name.endswith(("_mdr", "_dtr")):        # MDR/DTR-filtered year files belong to mdr_full_pipeline.py
                continue
            out.add(name[:-len("_baseline")] if name.endswith("_baseline") else name)
        return out
    return sorted(labels(a) & labels(b))


def year_file(year_dir, species, g, year):
    """Path of a year export: <age group>_<year>.xlsx, or the unfiltered <age group>_<year>_baseline.xlsx."""
    d = os.path.join(year_dir, safe(species))
    for suf in ("", "_baseline"):
        p = os.path.join(d, f"atlas_pair_{g}_{year}{suf}.xlsx")
        if os.path.exists(p):
            return p
    return os.path.join(d, f"atlas_pair_{g}_{year}.xlsx")


def gather_year_entries(species, drug, a, b, year_dir):
    """[(year, entry_a, entry_b)] for every year with >= MIN_GLOBAL_N isolates of this drug in both age groups."""
    entries = []
    for year in discover_years(species, a, b, year_dir):
        ea = parse_flat_absolute(year_file(year_dir, species, a, year)).get(drug)
        eb = parse_flat_absolute(year_file(year_dir, species, b, year)).get(drug)
        if ea and eb and ea["N"] >= MIN_GLOBAL_N and eb["N"] >= MIN_GLOBAL_N:
            entries.append((year, ea, eb))
    return entries


def needed_files(kind, species, a, b, comp_a, comp_b, site_folders, year_dir):
    """Files a flagged finding needs for its specimen source (kind="source") or year (kind="year") check.
    Levels = specimen sources (or years) with >= MIN_GLOBAL_N isolates of the species in BOTH age groups, read from
    the main exports. Returns [(age group, level, path to save to, already there?)]."""
    section = "Specimen Source Group" if kind == "source" else "Year"
    ca, cb = pooled_category_counts(comp_a, section), pooled_category_counts(comp_b, section)
    levels = sorted(str(lv) for lv in set(ca) & set(cb) if ca[lv] >= MIN_GLOBAL_N and cb[lv] >= MIN_GLOBAL_N)
    out = []
    for lv in levels:
        for g in (a, b):
            if kind == "source":
                name = f"atlas_pair_{g}_{safe_site(lv)}"
                folders = site_folders or [os.path.join(HERE, "data", "specimen_source")]
                have = any(os.path.exists(os.path.join(d, safe(species), name + suf + ".xlsx"))
                           for d in folders for suf in ("", "_baseline"))
                path = os.path.join(folders[0], safe(species), name + ".xlsx")
            else:
                path = year_file(year_dir, species, g, lv)
                have = os.path.exists(path)
            out.append((g, lv, path, have))
    return out


def files_missing_drug(drug, paths):
    """Existing exports in `paths` that do not contain `drug` at all (ATLAS can silently leave a drug out of an
    export when other filters are set). These need to be downloaded again with the drug ticked."""
    return [p for p in paths if os.path.exists(p) and drug not in parse_flat_absolute(p)]


def stratified_check(entries, unit_name):
    """Factor-adjusted shift: log2 shift within each level of the factor (specimen source or year), weighted by the
    harmonic mean of the two age groups' isolate counts, with a bootstrap that resamples the levels (2,000 draws).
    Same math as the country adjustment, with the factor's levels in place of countries.
    Returns (columns, survives) where survives is True, False, or "unknown" (fewer than 3 usable levels)."""
    cols = {"n_used": len(entries), "used": "; ".join(u for u, _, _ in entries)}
    if len(entries) < MIN_QUALIFYING_COUNTRIES:
        cols["check"] = (f"flagged, not checked (only {len(entries)} {unit_name}(s) with >= {MIN_GLOBAL_N} "
                         f"isolates in both age groups)")
        return cols, "unknown"
    point = weighted_site_shift(entries)
    boot = bootstrap_site_ci(entries)
    survives = bool(boot and boot["excludes_zero"])
    parts = []
    for u, ea, eb in entries:
        ga, gb = geomean_from_counts(ea["counts"]), geomean_from_counts(eb["counts"])
        if ga and gb:
            parts.append(f"{u}: {math.log2(gb / ga):+.3f}")
    cols.update(adjusted_log2_shift=round(point, 3) if point is not None else "",
                ci_lo=(boot or {}).get("ci_lo"), ci_hi=(boot or {}).get("ci_hi"), per_level_shifts="; ".join(parts),
                check=f"survives {unit_name} adjustment" if survives
                else f"weakened by {unit_name} adjustment (CI crosses zero)")
    return cols, survives


def parse_mic_classes(path, sheet="Absolute"):
    """{drug: {mic: [isolate count pooled over countries, S/I/R class or None]}} from ATLAS's cell colours."""
    cat = {COLOR_S.upper(): "S", COLOR_I.upper(): "I", COLOR_R.upper(): "R"}
    def mic(v):
        s = str(v).strip()
        for tok in ("\u2264", "<=", ">", "<"):
            s = s.replace(tok, "")
        try:
            return float(s)
        except ValueError:
            return None
    ws = openpyxl.load_workbook(path, data_only=True)[sheet]
    out = {}
    header, inside = None, False
    for row in ws.iter_rows(values_only=False):
        f0 = row[0].value
        if isinstance(f0, str) and ("isolates at each MIC" in f0 or "frequency distribution" in f0):
            header, inside = None, False
            continue
        if f0 == "Antibacterial":
            header, inside = [c.value for c in row], True
            continue
        if not inside:
            continue
        if not isinstance(f0, str) or not isinstance(row[1].value, (int, float)):
            inside = False
            continue
        for j, cell in enumerate(row[2:], start=2):
            if not cell.value:
                continue
            fill = cell.fill.fgColor.rgb if cell.fill and cell.fill.patternType else None
            e = out.setdefault(f0, {}).setdefault(mic(header[j]), [0, None])
            e[0] += int(cell.value)
            c = cat.get(fill.upper() if isinstance(fill, str) else fill)
            if c:
                e[1] = c
    return out


def median_mic_class(d):
    """(median MIC, its S/I/R class) of a pooled {mic: [count, class]} distribution."""
    tot, acc = sum(c for c, _ in d.values()), 0
    for m in sorted(d):
        acc += d[m][0]
        if acc >= tot / 2:
            return m, d[m][1]
    return None, None


def specimen_source_mix(comp):
    """{source: count} pooled over every country in one export."""
    return pooled_category_counts(comp, "Specimen Source Group")


# ============================================================================
# Permutation negative control (optional)
# ============================================================================
_FAMILY, _NBOOT = [], N_BOOTSTRAP_GATE


def _init_worker(family, n_boot):
    global _FAMILY, _NBOOT
    _FAMILY, _NBOOT = family, n_boot


def _one_permutation(k):
    rng = np.random.default_rng([PERM_SEED, k])
    nulls = [make_null_entries(e, rng) for e in _FAMILY]
    seeds = rng.integers(0, 2**31 - 1, size=len(_FAMILY))
    shifts = [weighted_shift(n) for n in nulls]
    hits = []
    for i, s in enumerate(shifts):
        if s is not None and abs(s) >= EFFECT_FLOOR:
            boot = bootstrap_ci(nulls[i], n_boot=_NBOOT, seed=int(seeds[i]))
            hits.append(bool(boot and boot["p_value"] < FDR_Q))
        else:
            hits.append(False)
    return hits


def permutation_control(family_entries, n_perm, workers):
    """Null rate per comparison of |shift| >= 0.5 and bootstrap p < 0.05 over n_perm null datasets (BH not applied,
    so this is an upper bound on the false positive rate)."""
    counts = np.zeros(len(family_entries), int)
    with Pool(workers, initializer=_init_worker, initargs=(family_entries, N_BOOTSTRAP_GATE)) as pool:
        for done, hits in enumerate(pool.imap_unordered(_one_permutation, range(n_perm)), 1):
            counts += np.array(hits, dtype=int)
            if done % 25 == 0 or done == n_perm:
                print(f"    permutations: {done}/{n_perm}", flush=True)
    return counts / n_perm


# ============================================================================
# Main run
# ============================================================================
ALL_FIELDS = [
    "species", "antibacterial", "comparison", "status", "direction", "log2_shift", "fold_change",
    "n_qualifying_countries", "bootstrap_ci_lo", "bootstrap_ci_hi", "bootstrap_p_value", "bh_q_value",
    "pct_resistant_a", "pct_resistant_b", "n_breakpoint_a", "n_breakpoint_b", "has_clinical_breakpoint",
    "expected_resistant_eucast", "median_mic_a", "median_class_a", "median_mic_b", "median_class_b",
    "median_class_change",
    "mixedlm_converged", "mixedlm_coef_log2_shift", "mixedlm_p_value",
    "specimen_source_tvd_pp", "specimen_source_differs",
    "specimen_source_check", "site_adjusted_log2_shift", "site_adjusted_ci_lo", "site_adjusted_ci_hi",
    "n_sites_used", "sites_used", "per_site_shifts", "robust_after_specimen_source_check",
    "year_tvd_pp", "year_differs",
    "year_check", "year_adjusted_log2_shift", "year_adjusted_ci_lo", "year_adjusted_ci_hi",
    "n_years_used", "years_used", "per_year_shifts", "robust_after_year_check",
    "robust_final", "permutation_null_rate",
]


def run(pairs=None, data_dir=DEFAULT_DATA_DIR, sample_size_table=DEFAULT_SAMPLE_SIZE_TABLE,
        site_folders=None, site_lists=None, year_dir=DEFAULT_YEAR_DIR,
        out_dir=DEFAULT_OUT_DIR, n_perm=0, workers=4):
    """pairs: {species: set of normalized drug names or None for all drugs}, or None for every species and drug."""
    t0 = time.time()
    if not os.path.isdir(data_dir):
        raise SystemExit(f"age_group_folder not found: {data_dir}\nSet age_group_folder in padi_input.txt to the folder "
                         f"that holds one subfolder per bacterium (see the instructions at the top of padi_input.txt).")
    if sample_size_table and not os.path.exists(sample_size_table):
        raise SystemExit(f"sample_size_table not found: {sample_size_table}\nRun pediatric_sample_size_check.py first.")
    for d in (site_folders or []):
        if not os.path.isdir(d):
            print(f"NOTE: specimen source folder not found, ignored: {d}", flush=True)
    if pairs is None:
        pairs = {d.replace("_", " "): None for d in sorted(os.listdir(data_dir))
                 if os.path.isdir(os.path.join(data_dir, d))}
    species_list = sorted(pairs)
    viable = load_viable_pairs(sample_size_table) if sample_size_table else None
    print(f"PADI full pipeline: {len(species_list)} species, "
          f"sample size filter {'on' if viable is not None else 'off'}", flush=True)

    rows, entries_by_key, files = [], {}, {}
    for species in species_list:
        d = os.path.join(data_dir, safe(species))
        paths = {g: os.path.join(d, f"atlas_pair_{g}.xlsx") for g in AGE_GROUPS}
        paths = {g: p for g, p in paths.items() if os.path.exists(p)}
        if len(paths) < 2:
            print(f"  {species}: {len(paths)} age group file(s) found in {d}, skipped "
                  f"(check the spelling matches the folder name)", flush=True)
            continue
        blocks = {g: parse_by_country_blocks(p) for g, p in paths.items()}
        files[species] = {"paths": paths, "blocks": blocks}
        panel = sorted({dr for bc in blocks.values() for dd in bc.values() for dr in dd})
        wanted_drugs = pairs[species]
        if wanted_drugs is not None:
            found = [dr for dr in panel if norm_name(dr) in wanted_drugs
                     or norm_name(EXPORT_TO_TABLE_NAME.get(dr, dr)) in wanted_drugs]
            found_names = {norm_name(dr) for dr in found} | {norm_name(EXPORT_TO_TABLE_NAME.get(dr, dr)) for dr in found}
            missing = sorted(wanted_drugs - found_names)
            if missing:
                print(f"  {species}: not found in the exports: {missing}. Drugs available: {', '.join(panel)}", flush=True)
            panel = found
        if viable is not None:
            dropped = [dr for dr in panel if not is_viable(viable, species, dr)]
            if dropped and wanted_drugs is not None:
                print(f"  {species}: failed the pediatric sample size check, not analyzed: {', '.join(dropped)}", flush=True)
            panel = [dr for dr in panel if is_viable(viable, species, dr)]
        n_tested = 0
        for drug in panel:
            for a, b in itertools.combinations(sorted(paths, key=AGE_GROUPS.index), 2):
                row = {f: "" for f in ALL_FIELDS} | {"species": species, "antibacterial": drug,
                                                    "comparison": f"{a}_vs_{b}"}
                entries = qualifying_countries(blocks[a], blocks[b], drug)
                row["n_qualifying_countries"] = len(entries)
                point = weighted_shift(entries) if entries else None
                if point is not None:
                    row.update(log2_shift=round(point, 3), fold_change=round(2 ** point, 3),
                               direction=f"higher in {b}" if point > 0 else f"higher in {a}")
                if len(entries) < MIN_QUALIFYING_COUNTRIES or point is None:
                    row["status"] = "insufficient_countries"
                    rows.append(row)
                    continue
                boot = bootstrap_ci(entries, n_boot=N_BOOTSTRAP_GATE)
                row.update(bootstrap_ci_lo=boot["ci_lo"], bootstrap_ci_hi=boot["ci_hi"],
                           bootstrap_p_value=boot["p_value"])
                entries_by_key[(species, drug, row["comparison"])] = entries
                rows.append(row)
                n_tested += 1
        print(f"  {species}: {len(panel)} drugs, {n_tested} comparisons bootstrapped "
              f"({(time.time() - t0) / 60:.1f} min)", flush=True)

    # ---- Benjamini-Hochberg across every bootstrapped comparison; robust = floor AND q < FDR_Q
    rows.sort(key=lambda r: (r["species"], r["antibacterial"], r["comparison"]))
    family = [r for r in rows if r["bootstrap_p_value"] != ""]
    qs = bh_qvalues([float(r["bootstrap_p_value"]) for r in family])
    for r, q in zip(family, qs):
        r["bh_q_value"] = q
        if abs(r["log2_shift"]) < EFFECT_FLOOR:
            r["status"] = "below_effect_floor"
        elif r["bootstrap_p_value"] >= FDR_Q:
            r["status"] = "not_significant"
        else:
            r["status"] = "robust" if q < FDR_Q else "fails_bh_fdr"
    robust = [r for r in family if r["status"] == "robust"]
    print(f"BH family m = {len(family)}; robust (|shift| >= {EFFECT_FLOOR}, q < {FDR_Q}): {len(robust)}", flush=True)

    # ---- Mixed model for every comparison that clears the floor (reported, not gating)
    floor_pass = [r for r in family if r["status"] != "below_effect_floor"]
    for r in floor_pass:
        a, b = r["comparison"].split("_vs_")
        mixed = fit_mixedlm(entries_by_key[(r["species"], r["antibacterial"], r["comparison"])], a, b)
        r.update(mixedlm_converged=mixed.get("converged"), mixedlm_coef_log2_shift=mixed.get("coef_log2_shift"),
                 mixedlm_p_value=mixed.get("p_value"))

    # ---- Follow-ups on candidates (robust + fails_bh_fdr)
    candidates = [r for r in family if r["status"] in ("robust", "fails_bh_fdr")]
    downloads = {}   # follow-up files a flagged finding still needs: (kind, species, age group, level) -> info
    site_folders = site_folders if site_folders is not None else [d for d, _ in DEFAULT_SITE_DIRS]
    site_dirs = site_search_order(site_folders)
    manifest = load_site_manifest(site_lists) if site_lists else None   # None: read sources from the file names
    comp_cache, bp_cache, cls_cache = {}, {}, {}
    for r in candidates:
        species, drug = r["species"], r["antibacterial"]
        a, b = r["comparison"].split("_vs_")
        f = files[species]
        r["expected_resistant_eucast"] = is_expected_resistant(species, drug)

        for g in (a, b):
            if (species, g) not in bp_cache:
                bp_cache[(species, g)] = parse_by_country_breakpoints(f["paths"][g])
                comp_cache[(species, g)] = parse_composition_sections(f["paths"][g])
        pr = country_adjusted_pct_resistant(f["blocks"][a], f["blocks"][b],
                                            bp_cache[(species, a)], bp_cache[(species, b)], drug)
        r["has_clinical_breakpoint"] = bool(pr)
        if pr:
            r.update(pct_resistant_a=pr[0], pct_resistant_b=pr[1], n_breakpoint_a=pr[2], n_breakpoint_b=pr[3])
            for g in (a, b):
                if (species, g) not in cls_cache:
                    cls_cache[(species, g)] = parse_mic_classes(f["paths"][g])
            ma, ca = median_mic_class(cls_cache[(species, a)].get(drug, {}))
            mb, cb = median_mic_class(cls_cache[(species, b)].get(drug, {}))
            r.update(median_mic_a=ma, median_class_a=ca, median_mic_b=mb, median_class_b=cb,
                     median_class_change=f"{ca} -> {cb}" if ca != cb else "no change")

        conf = composition_flags(comp_cache[(species, a)], comp_cache[(species, b)])
        spec = conf.get("Specimen Source Group", {})
        year = conf.get("Year", {})
        r.update(specimen_source_tvd_pp=spec.get("tvd_pp"), specimen_source_differs=spec.get("differs", False),
                 year_tvd_pp=year.get("tvd_pp"), year_differs=year.get("differs", False))

        # Specimen source check (only when flagged)
        if not spec.get("differs"):
            r.update(specimen_source_check="not needed (specimen source mix matches between age groups)",
                     robust_after_specimen_source_check=True)
        else:
            if manifest is not None:
                sites = sorted(set(manifest.get(species, {}).get(a, [])) & set(manifest.get(species, {}).get(b, [])))
            else:
                sites = discover_sources(species, a, b, site_folders)
            cols, ok = stratified_check(gather_site_entries(species, drug, a, b, sites, site_dirs), "specimen source")
            if ok == "unknown":
                need = needed_files("source", species, a, b, comp_cache[(species, a)], comp_cache[(species, b)],
                                    site_folders, year_dir)
                missing = [x for x in need if not x[3]]
                if missing:
                    for g, lv, path, _ in missing:
                        downloads.setdefault(("source", species, g, lv), {"path": path, "findings": set()})
                        downloads[("source", species, g, lv)]["findings"].add(f"{drug} {r['comparison']}")
                    cols["check"] = (f"flagged, files needed: download {len(missing)} specimen source file(s) "
                                     f"listed in files_to_download.txt, then run again")
                else:
                    existing = []
                    for g, lv, _, _ in need:
                        for d_, suf in site_dirs:
                            p = os.path.join(d_, safe(species), f"atlas_pair_{g}_{safe_site(lv)}{suf}.xlsx")
                            if os.path.exists(p):
                                existing.append((g, lv, p))
                                break
                    lacking = set(files_missing_drug(drug, [p for _, _, p in existing]))
                    for g, lv, p in existing:
                        if p in lacking:
                            key = ("source", species, g, lv)
                            downloads.setdefault(key, {"path": p, "findings": set(), "redownload": set()})
                            downloads[key].setdefault("redownload", set()).add(drug)
                            downloads[key]["findings"].add(f"{drug} {r['comparison']}")
                    if lacking:
                        cols["check"] = (f"flagged, files need re-downloading: {drug} is missing from "
                                         f"{len(lacking)} specimen source file(s) (see files_to_download.txt)")
            r.update(specimen_source_check=cols["check"], n_sites_used=cols["n_used"], sites_used=cols["used"],
                     site_adjusted_log2_shift=cols.get("adjusted_log2_shift", ""),
                     site_adjusted_ci_lo=cols.get("ci_lo", ""), site_adjusted_ci_hi=cols.get("ci_hi", ""),
                     per_site_shifts=cols.get("per_level_shifts", ""), robust_after_specimen_source_check=ok)

        # Year check (only when flagged)
        if not year.get("differs"):
            r.update(year_check="not needed (year mix matches between age groups)", robust_after_year_check=True)
        else:
            cols, ok = stratified_check(gather_year_entries(species, drug, a, b, year_dir), "year")
            if ok == "unknown":
                need = needed_files("year", species, a, b, comp_cache[(species, a)], comp_cache[(species, b)],
                                    site_folders, year_dir)
                missing = [x for x in need if not x[3]]
                if missing:
                    for g, lv, path, _ in missing:
                        downloads.setdefault(("year", species, g, lv), {"path": path, "findings": set()})
                        downloads[("year", species, g, lv)]["findings"].add(f"{drug} {r['comparison']}")
                    cols["check"] = (f"flagged, files needed: download {len(missing)} year file(s) "
                                     f"listed in files_to_download.txt, then run again")
            r.update(year_check=cols["check"], n_years_used=cols["n_used"], years_used=cols["used"],
                     year_adjusted_log2_shift=cols.get("adjusted_log2_shift", ""),
                     year_adjusted_ci_lo=cols.get("ci_lo", ""), year_adjusted_ci_hi=cols.get("ci_hi", ""),
                     per_year_shifts=cols.get("per_level_shifts", ""), robust_after_year_check=ok)

        checks = (r["robust_after_specimen_source_check"], r["robust_after_year_check"])
        r["robust_final"] = (False if r["status"] != "robust"
                             else True if all(c is True for c in checks)
                             else False if any(c is False for c in checks) else "unknown")

    # ---- Optional permutation negative control over the whole BH family
    if n_perm:
        print(f"Permutation control: {n_perm} null datasets over {len(family)} comparisons, {workers} workers", flush=True)
        keys = [(r["species"], r["antibacterial"], r["comparison"]) for r in family]
        rates = permutation_control([entries_by_key[k] for k in keys], n_perm, workers)
        for r, rate in zip(family, rates):
            r["permutation_null_rate"] = round(float(rate), 4)

    # ---- Write outputs
    os.makedirs(out_dir, exist_ok=True)
    def write(name, data):
        with open(os.path.join(out_dir, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=ALL_FIELDS, extrasaction="ignore")
            w.writeheader()
            w.writerows(data)
    write("padi_all_comparisons.csv", rows)
    write("padi_age_trends.csv", robust)

    # Specimen source mix per age group, for species with a flagged candidate
    mix_rows = []
    for species in sorted({r["species"] for r in candidates if r["specimen_source_differs"]}):
        for g in sorted(files[species]["paths"], key=AGE_GROUPS.index):
            if (species, g) not in comp_cache:
                comp_cache[(species, g)] = parse_composition_sections(files[species]["paths"][g])
            mix = specimen_source_mix(comp_cache[(species, g)])
            n = sum(mix.values())
            for src_name, cnt in sorted(mix.items(), key=lambda x: -x[1]):
                mix_rows.append({"species": species, "age_group": g, "specimen_source": src_name,
                                 "n_isolates": cnt, "pct": round(100 * cnt / n, 1) if n else ""})
    with open(os.path.join(out_dir, "padi_specimen_source_mix.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["species", "age_group", "specimen_source", "n_isolates", "pct"])
        w.writeheader()
        w.writerows(mix_rows)

    # Calibration of the 15 pp specimen source threshold
    retested = [r for r in candidates if r["site_adjusted_log2_shift"] not in ("", None) and r["specimen_source_tvd_pp"]]
    per15 = [abs(r["log2_shift"] - r["site_adjusted_log2_shift"]) / r["specimen_source_tvd_pp"] * 15 for r in retested]
    calibration = {"n_retested_within_specimen_sources": len(retested),
                   "max_abs_pooled_minus_adjusted_log2_per_15pp": round(max(per15), 3) if per15 else None}

    # Median MIC class changes: final robust findings with a clinical breakpoint, not expected resistant
    final_bp = [r for r in robust if r["robust_final"] is True and r["has_clinical_breakpoint"]
                and not r["expected_resistant_eucast"]]
    changes = {}
    for r in final_bp:
        changes[r["median_class_change"]] = changes.get(r["median_class_change"], 0) + 1
    median_classes = {"n_final_robust": sum(1 for r in robust if r["robust_final"] is True),
                      "n_final_robust_no_breakpoint": sum(1 for r in robust if r["robust_final"] is True
                                                          and not r["has_clinical_breakpoint"]),
                      "n_with_breakpoint_not_expected_resistant": len(final_bp),
                      "n_class_changed": sum(1 for r in final_bp if r["median_class_change"] != "no change"),
                      "n_moved_to_resistant": sum(1 for r in final_bp if r["median_class_b"] == "R"
                                                  and r["median_class_a"] != "R"),
                      "transitions": changes}
    status_counts = {}
    for r in rows:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1
    def tally(col):
        out = {}
        for r in robust:
            out[str(r[col])] = out.get(str(r[col]), 0) + 1
        return out
    summary = {"n_species": len(files), "n_comparisons": len(rows), "bh_family_m": len(family),
               "status_counts": status_counts, "n_robust": len(robust),
               "flagged_specimen_source": sum(1 for r in robust if r["specimen_source_differs"]),
               "flagged_year": sum(1 for r in robust if r["year_differs"]),
               "robust_after_specimen_source_check": tally("robust_after_specimen_source_check"),
               "robust_after_year_check": tally("robust_after_year_check"),
               "n_robust_final": tally("robust_final").get("True", 0),
               "specimen_source_calibration": calibration, "median_mic_classes": median_classes,
               "permutations": n_perm,
               "minutes": round((time.time() - t0) / 60, 1)}
    with open(os.path.join(out_dir, "padi_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"Wrote padi_all_comparisons.csv, padi_age_trends.csv, padi_specimen_source_mix.csv, "
          f"padi_summary.json to {out_dir}")

    # ---- Follow-up files still needed
    dl_path = os.path.join(out_dir, "files_to_download.txt")
    if not downloads:
        if os.path.exists(dl_path):
            os.remove(dl_path)
        return summary
    waiting = [r for r in robust if r["robust_final"] == "unknown"]
    lines = [
        "FILES TO DOWNLOAD FOR THE SPECIMEN SOURCE / YEAR CHECKS",
        "",
        f"{len(waiting)} robust finding(s) have a specimen source or year mix that differs between the age groups.",
        "They need the files below before they can be confirmed. Until then their robust_final is 'unknown'.",
        "Download each file from ATLAS (MIC Distribution), save it with exactly the name and folder shown,",
        "then run padi_full_pipeline.py again. One file covers every drug, so download it once.",
        "Files marked RE-DOWNLOAD already exist but are missing the drug shown: ATLAS can leave a drug out of",
        "the export when other filters are set, so tick that drug explicitly and replace the file.",
        "",
        "How to download each file:",
        "  - Select the bacterium shown and ALL antibiotics, ALL years (except for year files).",
        "  - Set the age filter to the age group shown.",
        "  - Turn View By: Country OFF (worldwide).",
        "  - Specimen source files: select ONLY the specimen source shown.",
        "  - Year files: select ONLY the year shown.",
        "  - Export, and save the .xlsx under the name below. In names, spaces are _ and '/' is -.",
        "",
    ]
    for kind, title in (("source", "SPECIMEN SOURCE FILES"), ("year", "YEAR FILES")):
        items = sorted((k, v) for k, v in downloads.items() if k[0] == kind)
        if not items:
            continue
        lines.append(f"==== {title} ====")
        for species in sorted({k[1] for k, _ in items}):
            sp_items = [(k, v) for k, v in items if k[1] == species]
            finds = sorted(set().union(*[v["findings"] for _, v in sp_items]))
            lines.append(f"{species}  (needed for: {'; '.join(finds)})")
            for (_, _, g, lv), v in sp_items:
                label = "specimen source" if kind == "source" else "year"
                redo = v.get("redownload")
                note = f"   RE-DOWNLOAD: missing {', '.join(sorted(redo))}" if redo else ""
                lines.append(f"  age {g:6s} {label} {lv:28s} -> {v['path']}{note}")
            lines.append("")
    with open(dl_path, "w") as fh:
        fh.write("\n".join(lines))
    n_files = len(downloads)
    print(f"\n>>> ACTION NEEDED: {len(waiting)} robust finding(s) need a specimen source or year check that this run "
          f"could not do, because {n_files} file(s) are missing or need re-downloading.\n>>> The list of files, where to save them and how "
          f"to download them is in {dl_path}\n>>> Download them, then run padi_full_pipeline.py again.")
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
    ap.add_argument("input", nargs="?", default=os.path.join(HERE, "padi_input.txt"), help="input file (default: padi_input.txt)")
    args = ap.parse_args()
    s = read_settings(args.input, need_bacteria=True)
    print(f"Input file: {s['path']}")
    run(pairs=s["pairs"], data_dir=require_folder(s, "age_group_folder", "age group"),
        sample_size_table=require_file(s, "sample_size_table", "pediatric_sample_size_check.py"),
        site_folders=s["specimen_source_folders"] or None, site_lists=[],
        year_dir=s["year_folder"] or os.path.join(os.path.dirname(s["path"]), "data", "year"),
        out_dir=s["results_folder"], n_perm=s["permutations"], workers=s["workers"])


if __name__ == "__main__":
    main()
