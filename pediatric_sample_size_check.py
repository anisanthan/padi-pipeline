"""
Pediatric sample size check (step 1). Reads the pediatric (0-18) export for every country and combines them into
one table over the FULL universe (ALL_SPECIES x ALL_DRUGS below =
266 x 48 = 12,768 pairs) and flags anything the passes disagree on.

For each (species, drug): per-country pediatric N is taken as the MAX across
passes (a pair one pass missed and another saw still counts), then the usual
screen applies: >=3 countries with N>=20 -> viable, else not_viable; pairs with
no rows anywhere -> no_rows.

Usage:  python3 pediatric_sample_size_check.py            (reads sample_size_input.txt: pediatric_folder, results_folder)
        python3 pediatric_sample_size_check.py my_input.txt
The country files go in pediatric_folder as atlas_ped_0to18_<Country>.xlsx (see sample_size_input.txt for how to download
them). If you downloaded in several rounds, put each round in a subfolder named pass1, pass2, ... instead.

Outputs:
  <results_folder>/pediatric_sample_size.csv        one row per pair: viable / not_viable / no_rows
  <results_folder>/pass_differences.csv             only if several download rounds disagree
Also prints species/drugs that appear in exports but not in the master lists
(spelling aliases or a truncated list) and countries missing a pass.
"""
import argparse
import csv
import glob
import os
import re
from collections import defaultdict

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
SUF = re.compile(r"\s*\(N = \d+\)\s*$")
PREFIX = "atlas_ped_0to18_"
# ATLAS spells some things differently from the master lists: drugs with hyphens ("Piperacillin-tazobactam")
# via UI_TO_OURS, and this species typo. Without this those pairs would silently show up as no_rows.
SPECIES_ALIAS = {"Enterobacter hormaechi": "Enterobacter hormaechei"}

MIN_COUNTRY_N = 20              # pediatric isolates per country
MIN_QUALIFYING_COUNTRIES = 3    # countries needed

# The pair universe: every species and every drug in ATLAS's pediatric (0-18) MIC panel, copied by hand from the
# ATLAS panel on 2026-09-19 (266 species, 48 drugs).
ALL_SPECIES = [
    'Acidaminococcus fermentans', 'Acinetobacter anitratus',
    'Acinetobacter baumannii', 'Acinetobacter baylyi',
    'Acinetobacter beijerinckii', 'Acinetobacter bereziniae',
    'Acinetobacter calcoaceticus', 'Acinetobacter colistiniresistens',
    'Acinetobacter courvalinii', 'Acinetobacter dijkshoorniae',
    'Acinetobacter dispersus', 'Acinetobacter guillouiae',
    'Acinetobacter gyllenbergii', 'Acinetobacter haemolyticus',
    'Acinetobacter indicus', 'Acinetobacter johnsonii',
    'Acinetobacter junii', 'Acinetobacter lactucae',
    'Acinetobacter lwoffii', 'Acinetobacter nosocomialis',
    'Acinetobacter parvus', 'Acinetobacter pitii',
    'Acinetobacter radioresistens', 'Acinetobacter schindleri',
    'Acinetobacter seifertii', 'Acinetobacter soli',
    'Acinetobacter spp', 'Acinetobacter tjernbergiae',
    'Acinetobacter towneri', 'Acinetobacter ursingii',
    'Acinetobacter variabilis', 'Acinetobacter, non-speciated',
    'Anaerococcus hydrogenalis', 'Anaerococcus murdochii',
    'Anaerococcus prevotii', 'Anaerococcus spp',
    'Anaerococcus vaginalis', 'Bacillus spp',
    'Bacteroides caccae', 'Bacteroides capsillosis',
    'Bacteroides eggerthii', 'Bacteroides faecis',
    'Bacteroides fragilis', 'Bacteroides merdeae',
    'Bacteroides nordii', 'Bacteroides ovatus',
    'Bacteroides salersyae', 'Bacteroides spp',
    'Bacteroides stercosis', 'Bacteroides thetaiotaomicron',
    'Bacteroides uniformis', 'Bacteroides vulgatus',
    'Bordetella trematum', 'Burkholderia cenocepacia',
    'Burkholderia cepacia', 'Citrobacter amalonaticus',
    'Citrobacter braakii', 'Citrobacter diversus',
    'Citrobacter farmeri', 'Citrobacter freundii',
    'Citrobacter gillenii', 'Citrobacter koseri',
    'Citrobacter sedlakii', 'Citrobacter spp',
    'Citrobacter youngae', 'Citrobacter, non-speciated',
    'Clostridioides difficile', 'Clostridioides spp',
    'Clostridium aldenense', 'Clostridium barati',
    'Clostridium bolteae', 'Clostridium butyricum',
    'Clostridium cadaveris', 'Clostridium clostridiiformis',
    'Clostridium colicanis', 'Clostridium hathewayi',
    'Clostridium histolyticum', 'Clostridium innocuum',
    'Clostridium limosum', 'Clostridium paraputrificum',
    'Clostridium perfringens', 'Clostridium ramosum',
    'Clostridium scindens', 'Clostridium septicum',
    'Clostridium sordellii', 'Clostridium sphenoides',
    'Clostridium sporogenes', 'Clostridium spp',
    'Clostridium subterminale', 'Clostridium symbiosum',
    'Clostridium tertium', 'Comamonas kerstersii',
    'Cronobacter sakazakii', 'Cutibacterium acnes',
    'Cutibacterium spp', 'Enterobacter agglomerans',
    'Enterobacter asburiae', 'Enterobacter bugandensis',
    'Enterobacter cancerogenus', 'Enterobacter cloacae',
    'Enterobacter hormaechei', 'Enterobacter intermedium',
    'Enterobacter kobei', 'Enterobacter ludwigii',
    'Enterobacter roggenkampii', 'Enterobacter sakazakii',
    'Enterobacter spp', 'Enterobacter xiangfangensis',
    'Enterobacter, non-speciated', 'Enterococcus avium',
    'Enterococcus casseliflavus', 'Enterococcus durans',
    'Enterococcus faecalis', 'Enterococcus faecium',
    'Enterococcus gallinarum', 'Enterococcus Group D',
    'Enterococcus hirae', 'Enterococcus raffinosus',
    'Enterococcus spp', 'Enterococcus, non-speciated',
    'Escherichia coli', 'Escherichia hermanii',
    'Escherichia spp', 'Eubacterium aerofaciens',
    'Eubacterium limosum', 'Finegoldia magna',
    'Finegoldia spp', 'Haemophilus influenzae',
    'Haemophilus parahaemolyticus', 'Haemophilus parainfluenzae',
    'Haemophilus spp', 'Hungatella hathewayi',
    'Klebsiella aerogenes', 'Klebsiella ornithinolytica',
    'Klebsiella oxytoca', 'Klebsiella ozaenae',
    'Klebsiella planticola', 'Klebsiella pneumoniae',
    'Klebsiella spp', 'Klebsiella variicola',
    'Klebsiella, non-speciated', 'Lelliottia amnigena',
    'Moraxella catarrhalis', 'Moraxella spp',
    'Morganella morganii', 'Morganella spp',
    'Neisseria gonorrhoeae', 'Pantoea agglomerans',
    'Pantoea dispersa', 'Pantoea spp',
    'Parabacteroides distasonis', 'Parabacteroides goldsteinii',
    'Parvimonas micra', 'Parvimonas spp',
    'Peptoniphilus asaccharolyticus', 'Peptoniphilus coxii',
    'Peptoniphilus grossensis', 'Peptoniphilus harei',
    'Peptoniphilus indolicus', 'Peptoniphilus lacrimalis',
    'Peptoniphilus spp', 'Peptostreptococcus anaerobius',
    'Peptostreptococcus hydrogenalis', 'Peptostreptococcus lactolyticus',
    'Peptostreptococcus magnus', 'Peptostreptococcus spp',
    'Phocaeicola vulgatus', 'Pluralibacter gergoviae',
    'Prevotella bergensis', 'Prevotella bivia',
    'Prevotella buccae', 'Prevotella buccalis',
    'Prevotella corporis', 'Prevotella denticola',
    'Prevotella disiens', 'Prevotella heparinolytica',
    'Prevotella histicola', 'Prevotella intermedia',
    'Prevotella loescheii', 'Prevotella maculosa',
    'Prevotella melaninogenica', 'Prevotella nanceinsis',
    'Prevotella nigrescens', 'Prevotella oralis',
    'Prevotella oris', 'Prevotella oulorum',
    'Prevotella salivae', 'Prevotella spp',
    'Prevotella timonensis', 'Proteus hauseri',
    'Proteus mirabilis', 'Proteus penneri',
    'Proteus spp', 'Proteus vulgaris',
    'Providencia alcalifaciens', 'Providencia rettgeri',
    'Providencia spp', 'Providencia stuartii',
    'Providencia, non-speciated', 'Pseudomonas aeruginosa',
    'Pseudomonas guariconensis', 'Pseudomonas monteilii',
    'Pseudomonas mosselii', 'Pseudomonas nitroreducens',
    'Pseudomonas otitidis', 'Pseudomonas putida',
    'Pseudomonas spp', 'Pseudomonas stutzeri',
    'Raoultella ornithinolytica', 'Raoultella planticola',
    'Salmonella spp', 'Serratia ficaria',
    'Serratia fonticola', 'Serratia liquefaciens',
    'Serratia marcescens', 'Serratia odorifera',
    'Serratia plymuthica', 'Serratia rubidaea',
    'Serratia spp', 'Serratia ureilytica',
    'Serratia, non-speciated', 'Staphylococcus argenteus',
    'Staphylococcus arlettae', 'Staphylococcus aureus',
    'Staphylococcus auricularis', 'Staphylococcus capitis',
    'Staphylococcus caprae', 'Staphylococcus Coagulase Negative',
    'Staphylococcus cohnii', 'Staphylococcus epidermidis',
    'Staphylococcus haemolyticus', 'Staphylococcus hominis',
    'Staphylococcus lugdunensis', 'Staphylococcus pasteuri',
    'Staphylococcus petrasii', 'Staphylococcus pettenkoferi',
    'Staphylococcus pseudointermedius', 'Staphylococcus saccharolyticus',
    'Staphylococcus saprophyticus', 'Staphylococcus schleiferi',
    'Staphylococcus sciuri', 'Staphylococcus simulans',
    'Staphylococcus spp', 'Staphylococcus warneri',
    'Staphylococcus xylosus', 'Stenotrophomonas maltophilia',
    'Streptococcus agalactiae', 'Streptococcus anginosus',
    'Streptococcus bovis', 'Streptococcus canis',
    'Streptococcus castoreus', 'Streptococcus constellatus',
    'Streptococcus dysgalactiae', 'Streptococcus gallolyticus',
    'Streptococcus intermedius', 'Streptococcus lutetiensis',
    'Streptococcus mitis_oralis', 'Streptococcus oralis',
    'Streptococcus parasanguinis', 'Streptococcus peroris',
    'Streptococcus pneumoniae', 'Streptococcus pyogenes',
    'Streptococcus salivarius', 'Streptococcus spp',
    'Streptococcus, Beta Hemolytic', 'Streptococcus, viridans group',
]

ALL_DRUGS = [
    'Amikacin', 'Amoxy/clav', 'Ampicillin', 'Ampicillin sulbactam',
    'Azithromycin', 'Aztreonam', 'Aztreonam avibactam', 'Cefepime',
    'Cefiderocol', 'Cefixime', 'Cefoperazone sulbactam', 'Cefoxitin',
    'Cefpodoxime', 'Ceftaroline', 'Ceftazidime', 'Ceftazidime avibactam',
    'Ceftibuten', 'Ceftolozane tazobactam', 'Ceftriaxone', 'Ciprofloxacin',
    'Clarithromycin', 'Clindamycin', 'Colistin', 'Daptomycin',
    'Doripenem', 'Ertapenem', 'Erythromycin', 'Gatifloxacin',
    'Gentamicin', 'Imipenem', 'Levofloxacin', 'Linezolid',
    'Meropenem', 'Meropenem vaborbactam', 'Metronidazole', 'Minocycline',
    'Moxifloxacin', 'Oxacillin', 'Penicillin', 'Pip/taz',
    'Quinupristin dalfopristin', 'Sulbactam', 'Teicoplanin', 'Tetracycline',
    'Tigecycline', 'Tobramycin', 'Trimethoprim sulfa', 'Vancomycin',
]

# ATLAS export spelling -> spelling in ALL_DRUGS (only where they differ)
UI_TO_OURS = {
    'Trimethoprim-sulfamethoxazole': 'Trimethoprim sulfa',
    'Amoxicillin-clavulanate': 'Amoxy/clav',
    'Piperacillin-tazobactam': 'Pip/taz',
    'Ampicillin-sulbactam': 'Ampicillin sulbactam',
    'Cefoperazone-sulbactam': 'Cefoperazone sulbactam',
    'Ceftazidime-avibactam': 'Ceftazidime avibactam',
    'Ceftolozane-tazobactam': 'Ceftolozane tazobactam',
    'Aztreonam-avibactam': 'Aztreonam avibactam',
    'Meropenem-vaborbactam': 'Meropenem vaborbactam',
    'Quinupristin-dalfopristin': 'Quinupristin dalfopristin',
}


def parse_species_blocks(path, sheet="Absolute"):
    """One export, all species stacked. Returns {species: {drug: N}} --
    just the isolate count per drug, not full MIC counts (all we need here)."""
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(min_row=1, max_row=ws.max_row, values_only=True))
    data, current, in_table = {}, None, False
    for row in rows:
        first = row[0]
        if first and isinstance(first, str) and "frequency distribution for" in first:
            current = first.split("frequency distribution for")[1].strip()
            data.setdefault(current, {})
            in_table = False
            continue
        if first and isinstance(first, str) and "isolates at each MIC" in first and " for " in first:
            current = first.rsplit(" for ", 1)[1].strip()
            data.setdefault(current, {})
            in_table = False
            continue
        if first == "Antibacterial":
            in_table = True
            continue
        if not in_table or current is None:
            continue
        if first is None:
            in_table = False
            continue
        if not isinstance(first, str) or not isinstance(row[1], (int, float)):
            in_table = False
            continue
        data[current][first] = int(row[1])
    return data


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
    ap.add_argument("input", nargs="?", default=os.path.join(HERE, "sample_size_input.txt"), help="input file (default: sample_size_input.txt)")
    args = ap.parse_args()
    s = read_settings(args.input)
    ped = require_folder(s, "pediatric_folder", "pediatric (0-18) country")
    out = s["results_folder"]
    os.makedirs(out, exist_ok=True)
    pass_dirs = sorted(d for d in glob.glob(os.path.join(ped, "pass*")) if os.path.isdir(d)) or [ped]
    if not any(glob.glob(os.path.join(d, PREFIX + "*.xlsx")) for d in pass_dirs):
        raise SystemExit(f"No {PREFIX}<Country>.xlsx files found in {ped}. See sample_size_input.txt for how to download "
                         f"and name them.")
    per_pass = {}                      # pass -> {(species, drug): {country: N}}
    countries_by_pass = {}
    for d in pass_dirs:
        counts, countries = defaultdict(dict), set()
        for path in sorted(glob.glob(os.path.join(d, PREFIX + "*.xlsx"))):
            country = os.path.basename(path)[len(PREFIX):-len(".xlsx")].replace("_", " ")
            countries.add(country)
            for sp, drugs in parse_species_blocks(path).items():
                sp = SUF.sub("", sp).strip()
                sp = SPECIES_ALIAS.get(sp, sp)
                for drug, n in drugs.items():
                    counts[(sp, UI_TO_OURS.get(drug, drug))][country] = n
        per_pass[os.path.basename(d)] = counts
        countries_by_pass[os.path.basename(d)] = countries

    all_countries = set().union(*countries_by_pass.values())
    for p, cs in countries_by_pass.items():
        gap = sorted(all_countries - cs)
        print(f"{p}: {len(cs)} countries" + (f" | MISSING vs other passes: {gap}" if gap else ""))

    merged = defaultdict(dict)
    diffs = []
    keys = set().union(*[set(c) for c in per_pass.values()])
    for k in keys:
        for country in set().union(*[set(c.get(k, {})) for c in per_pass.values()]):
            vals = {p: c.get(k, {}).get(country) for p, c in per_pass.items()}
            present = [v for v in vals.values() if v is not None]
            merged[k][country] = max(present)
            if len(present) < len(per_pass) or len(set(present)) > 1:
                diffs.append({"species": k[0], "antibacterial": k[1], "country": country,
                              **{p: ("" if v is None else v) for p, v in vals.items()}})

    universe = [(s, d) for s in ALL_SPECIES for d in ALL_DRUGS]
    # Also classify every pair ATLAS actually returned that the master lists don't cover (e.g. Tebipenem,
    # Ceftaroline-avibactam, Ceftibuten-avibactam, Staphylococcus pseudointermedius), so nothing viable is dropped.
    master_set = set(universe)
    universe += sorted(k for k in keys if k not in master_set)
    rows = []
    for k in universe:
        by_c = merged.get(k, {})
        q = sum(1 for n in by_c.values() if n >= MIN_COUNTRY_N)
        status = "no_rows" if not by_c else ("viable" if q >= MIN_QUALIFYING_COUNTRIES else "not_viable")
        rows.append({"species": k[0], "antibacterial": k[1], "in_master_universe": k in master_set, "status": status,
                     "n_qualifying_countries": q, "n_countries_with_rows": len(by_c),
                     "total_pediatric_n": sum(by_c.values())})
    with open(os.path.join(out, "pediatric_sample_size.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    if diffs:
        with open(os.path.join(out, "pass_differences.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(diffs[0])); w.writeheader(); w.writerows(diffs)

    from collections import Counter
    c = Counter(r["status"] for r in rows)
    print(f"\n{len(universe)} pairs classified ({len(master_set)} master + {len(universe) - len(master_set)} extra from exports): "
          + ", ".join(f"{k}={v}" for k, v in c.most_common()))
    print("viable, master pairs:", sum(1 for r in rows if r["status"] == "viable" and r["in_master_universe"]),
          "| viable, extra pairs:", sum(1 for r in rows if r["status"] == "viable" and not r["in_master_universe"]))
    print(f"passes disagreeing on {len(diffs)} (pair, country) cells"
          + (f" -> {os.path.join(out, 'pass_differences.csv')}" if diffs else " (none)"))
    seen_sp, seen_dr = {s for s, _ in keys}, {d for _, d in keys}
    print("species in exports but not in master:", sorted(seen_sp - set(ALL_SPECIES)))
    print("drugs in exports but not in master:  ", sorted(seen_dr - set(ALL_DRUGS)))
    print(f"wrote {os.path.join(out, 'pediatric_sample_size.csv')}")


if __name__ == "__main__":
    main()
