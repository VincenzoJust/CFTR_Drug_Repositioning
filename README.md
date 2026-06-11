# QSAR Pipeline for Drug Repositioning in Cystic Fibrosis (CFTR)

Computational drug-repositioning study targeting the **CFTR** protein for cystic
fibrosis. A QSAR Random Forest model (trained on a fused ChEMBL + Papyrus +
BindingDB bioactivity dataset) is combined with molecular-docking energies into a
single integrated score that ranks repurposable, already-approved drugs.

> Master's Thesis (TFM) — MSc in Bioinformatics and Biostatistics (UOC, 2025–2026)
> Author: **Vicente Valero Just**

---

## What this notebook does

This repository contains the **QSAR + integration** block of the project. The
full pipeline has three blocks; this notebook covers blocks 2–3 and consumes the
output of block 1:

1. **PPI network analysis** (STRING + DrugBank) — identifies hub proteins
   (CFTR, HSP90AA1, PRKACA, HSPA8) and the 7 repositioning candidates. *Done
   externally; its output is the candidate list hard-coded in section 9.*
2. **Molecular docking** (AutoDock Vina) on two CFTR structures: PDB **6MSM**
   (apo WT, cryo-EM) and an **AlphaFold2 ΔF508** model. *Done externally; the
   resulting binding energies are hard-coded in section 11.*
3. **QSAR model** (this notebook): a Random Forest trained on the fused dataset
   predicts pIC50 for each molecule. Docking and QSAR are then merged into an
   **integrated score = 0.60 × docking + 0.40 × QSAR**.

The model is validated with positive controls (Ivacaftor and the approved CFTR
correctors) and negative controls (drugs with no known CFTR activity), and each
prediction is flagged with an **applicability-domain** label based on Tanimoto
similarity to the training set.

---

## Repository structure

```
.
├── CFTR_QSAR_pipeline.ipynb   # main notebook (run top to bottom)
├── requirements.txt
├── data/
   ├── README.md              # how to obtain the input data
   └── bindingdb_cftr.tsv     # you must download this (see data/README.md)

```

---

## Input data

| Source | How it is obtained | Manual download? |
|--------|--------------------|------------------|
| **ChEMBL** (target CHEMBL4051) | Downloaded automatically via the ChEMBL API | No |
| **Papyrus** v05.7 (UniProt P13569) | Downloaded automatically and cached (~several GB on first run) | No |
| **BindingDB** (target CFTR) | TSV exported manually from the website | **Yes** → `data/bindingdb_cftr.tsv` |

Only **BindingDB** requires a manual step. Full instructions are in
[`data/README.md`](data/README.md). If the TSV is missing, the notebook still
runs using ChEMBL + Papyrus only and prints a warning.

The docking energies (section 11) and PLIP interaction counts (section 16) are
**hard-coded** results from the external docking/PLIP work, so the final ranking
can be reproduced without re-running AutoDock Vina.

---

## How to run

```bash
# 1. (recommended) create an environment
python -m venv .venv && source .venv/bin/activate

# 2. install dependencies
pip install -r requirements.txt

# 3. (optional) place data/bindingdb_cftr.tsv  — see data/README.md

# 4. open and run the notebook top to bottom
jupyter lab CFTR_QSAR_pipeline.ipynb
```

First run is slow because Papyrus downloads and caches a large file. Later runs
use the cache.

---

## Outputs (written to `results/`)

| File | Content |
|------|---------|
| `dataset_full.csv` | Merged, deduplicated training dataset (SMILES + pIC50 + source) |
| `model_comparison.csv` | ChEMBL-only vs fused-dataset model metrics |
| `final_ranking.csv` | Final candidate ranking with score, docking, pIC50, AD |
| `control_validation.csv` | All molecules with predicted pIC50, docking and AD |
| `fig_final_ranking.png` | Ranking bar chart colored by applicability domain |
| `fig_control_validation.png` | Predicted pIC50 vs docking, by group |
| `fig_feature_importance.png` | Top-15 Random Forest descriptors |
| `fig_plip_comparison.png` | PLIP interaction profile of the top-4 candidates |

---

## Method notes and caveats

- The fused dataset (ChEMBL + Papyrus + BindingDB) is more stable than ChEMBL
  alone; the near-equal R² despite more data suggests the public CFTR
  bioactivity space is close to its informational limit with conventional
  descriptors.
- All candidates fall in the **marginal applicability domain** (Tanimoto
  ≈ 0.23–0.33), which is why docking is weighted higher (60%) than QSAR (40%).
- The integrated score is **relative to the candidate group** (MinMax-normalized):
  0.000 means "worst of this set", not "no affinity".
- Docking score differences of ~0.3 kcal/mol fall within Vina's error range
  (~2–3 kcal/mol) and should not be read as one drug "beating" another.
- Fostamatinib is a **prodrug**; its docking/PLIP results are interpreted with
  caution and its active metabolite **R406** is included as a check.

---

## Reference structures

- **PDB 6MSM** — apo wild-type CFTR (cryo-EM); the docking receptor.
- **PDB 6O2P** — CFTR in complex with Ivacaftor (used as a structural reference).
- **AlphaFold2 ΔF508** — model of the F508del mutant.

---

## Tools

RDKit · scikit-learn · AutoDock Vina · PrankWeb · AlphaFold2 · PyMOL · PLIP ·
SwissADME · STRING · DrugBank

## License

Released under the MIT License (see `LICENSE`). Bioactivity data remain subject
to the original ChEMBL, Papyrus and BindingDB licenses.
