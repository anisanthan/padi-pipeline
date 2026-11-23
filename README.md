# PADI: Pediatric Adult Divergence of Antimicrobial Resistance

Code for the paper "Pediatric Adult Divergence (PADI) of Antimicrobial Resistance (AMR) across Bacteria in ATLAS".
It compares minimum inhibitory concentrations (MICs) and multidrug resistance between children (0-18), adults (19-64)
and older adults (65+) using data from Pfizer's ATLAS surveillance database.

The ATLAS data is not included. Each script reads Excel files you download from
[ATLAS Interactive](https://atlas-surveillance.com) (free, registration required). The input file for each step
explains exactly which files to download, which filters to use, and how to name them.

## Requirements

Python 3.10 or later with numpy, scipy, pandas, openpyxl and statsmodels:

```
pip install numpy scipy pandas openpyxl statsmodels
```

## How to run

Run the five scripts in this order. Each one has its own input file: open it, follow the instructions at the top,
save it, then run the script from this folder.

| Step | Script | Input file | What it does |
|---|---|---|---|
| 1 | `pediatric_sample_size_check.py` | `sample_size_input.txt` | Finds the bacteria-drug pairs with at least 20 pediatric isolates in 3 or more countries. |
| 2 | `padi_full_pipeline.py` | `padi_input.txt` | The PADI method: country-adjusted MIC shift between age groups, bootstrap, Benjamini-Hochberg correction, and specimen source and year checks. |
| 3 | `mdr_full_pipeline.py` | `mdr_input.txt` | Compares the share of MDR and DTR isolates between age groups, with the same country adjustment and checks. |
| 4 | `marker_pipeline.py` | `marker_input.txt` | Tests whether MRSA and CRAB explain the S. aureus and A. baumannii age shifts. Needs step 2. |
| 5 | `mdr_mic_correlation.py` | `correlation_input.txt` | Correlates each bacterium's MDR gap with its median MIC shift. Needs steps 2 and 3. |

```
python3 pediatric_sample_size_check.py
python3 padi_full_pipeline.py
python3 mdr_full_pipeline.py
python3 marker_pipeline.py
python3 mdr_mic_correlation.py
```

To analyze only some bacteria or drugs, list them at the bottom of `padi_input.txt` (and `mdr_input.txt`).
The Benjamini-Hochberg correction covers every comparison in a run, so the paper's results come from analyzing
all bacteria and drugs (`all`), and a shorter list gives different q-values.

## Specimen source and year checks

You do not need specimen source or year files before the first run. If a finding's specimen source or year mix
differs between age groups, steps 2 and 3 end with **ACTION NEEDED** and write `files_to_download.txt` in their
results folder, listing exactly which files to download and where to save them. Download them and run the step again.

## Folder layout

By default the input files expect this layout (change the paths in the input files to use another):

```
data/
  pediatric_by_country/    step 1: atlas_ped_0to18_<Country>.xlsx
  age_groups/              steps 2 and 4: <Bacterium>/atlas_pair_<age group>.xlsx
  mdr_dtr/                 step 3: <Bacterium>/atlas_pair_<age group>_<baseline|mdr|dtr>.xlsx
  specimen_source/         steps 2 and 3, only if asked for
  year/                    steps 2 and 3, only if asked for
results/                   created by the scripts, one subfolder per step
```

## Citation

If you use this code, please cite the paper (citation to be added).
