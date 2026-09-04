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
3. **QSAR model** (this notebook): an ensemble (Random Forest + XGBoost + SVR
   with a Tanimoto kernel) trained on the fused dataset predicts activity for
   each molecule, with a per-prediction uncertainty interval. Docking and QSAR
   are then merged into an **integrated score = 0.60 × docking + 0.40 × QSAR**.

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
├── docs/
│   └── IMPROVEMENTS.md        # what was changed, why, and the measured effect
├── data/
│   ├── README.md              # how to obtain the input data
│   ├── bindingdb_cftr.tsv     # you must download this (see data/README.md)
│   └── cache/                 # frozen source downloads (committed, ~0.5 MB)
└── results/                   # generated CSVs
    └── figures/               # generated figures (PNG, 300 dpi)
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

### A fresh clone needs no large download

`data/cache/` contains frozen copies of exactly the data this notebook consumes —
the ChEMBL pull (208 KB) and the CFTR high-quality slice of Papyrus (248 KB) —
and both are committed. The notebook is **cache-first**, so a clean checkout runs
end to end offline, in seconds, without touching the ChEMBL API or downloading
the 1.3 GB Papyrus dataset. This is verified by running with `~/.data/papyrus`
removed entirely.

To re-download from source instead, set `REFRESH_DATA = True` in the setup cell.
That refreshes both caches and their provenance sidecars. It is also required if
you want to change the Papyrus quality threshold, since the cache is stored
*after* the `Quality == 'High'` filter.

Keeping the data frozen is deliberate: ChEMBL and Papyrus both change over time,
and a thesis result should not shift between runs.

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

First run is slow because Papyrus downloads and caches a large file (~10 GB in
`~/.data/papyrus`). Later runs use the cache.

> **If you use conda:** some conda-forge MKL builds (seen with `mkl 2026.1.0`)
> crash the Python interpreter on *any* numpy matrix multiply, which silently
> breaks the neural-network comparison. If `numpy.random.rand(200,200) @ .T`
> kills your Python, switch the BLAS:
> `conda install -c conda-forge "libblas=*=*openblas"`

---

## Outputs (written to `results/`)

| File | Content |
|------|---------|
| `dataset_full.csv` | Merged, deduplicated training dataset (SMILES, pActivity, replicate spread, endpoint, source) |
| `model_comparison.csv` | ChEMBL-only vs fused-dataset model metrics |
| `model_split_comparison.csv` | Every model under random / scaffold / Butina cross-validation |
| `final_ranking.csv` | Candidate ranking with score, docking, predicted activity + 90% interval, AD |
| `control_validation.csv` | All molecules with prediction, interval, AD, and whether they are in the training set |
| `data_provenance.json` | Which database versions produced these results (ChEMBL release, Papyrus version, BindingDB SHA256) |
| `figures/fig_final_ranking.png` | Ranking bar chart + the prediction intervals behind it |
| `figures/fig_split_comparison.png` | R² per model under each split — the generalisation gap |
| `figures/fig_control_validation.png` | Predicted activity vs docking, by group |
| `figures/fig_feature_importance.png` | Top-15 Random Forest features |
| `figures/fig_plip_comparison.png` | PLIP interaction profile of the top-4 candidates |

---

## Method notes and caveats

**Reported performance.** The headline figure is the **scaffold-split R² = 0.640**
(Bemis-Murcko `GroupKFold`), not the random 5-fold R² of 0.758. Bioactivity
datasets are dominated by analogue series, so a random split tests the model on
molecules nearly identical to ones it has memorised. Since all seven candidates
are structurally unrelated to the training data (max Tanimoto 0.26–0.33), the
scaffold split is the regime that *matches* how the model is being used. A
y-scrambling control returns R² = −0.10, confirming the signal is real.

**The QSAR term cannot rank the candidates.** Their predicted activities span
0.374 units while the uncertainty on any single prediction is 0.22–0.56. The
ordering is therefore driven almost entirely by docking — which the weight
sensitivity analysis independently confirms. The 60/40 weights were deliberately
*not* re-tuned after seeing this; instead every prediction is reported with a
90% interval and the limitation is stated explicitly.

**Controls are read by mechanism.** Ivacaftor is a potentiator; Lumacaftor,
Tezacaftor and Elexacaftor are correctors. About 87% of the training endpoint is
potentiator-style functional EC50, so correctors are genuinely weak *in this
assay class* (Lumacaftor's own measured median is pEC50 5.59). Only the
potentiator is required to clear the negative controls.

**Other caveats.**

- All candidates fall in the **marginal applicability domain** (Tanimoto
  ≈ 0.26–0.33), which is why docking is weighted higher (60%) than QSAR (40%).
- The integrated score is **relative to the candidate group** (MinMax-normalized):
  0.000 means "worst of this set", not "no affinity".
- Docking score differences of ~0.3 kcal/mol fall within Vina's error range
  (~2–3 kcal/mol) and should not be read as one drug "beating" another.
- Fostamatinib is a **prodrug**; its docking/PLIP results are interpreted with
  caution and its active metabolite **R406** is included as a check.
- The `pIC50 >= 5` activity filter is retained, so the model never sees an
  inactive compound and cannot predict below ~5.

A full account of what was changed, why, and the measured effect of each change
— including what was tried and did **not** work — is in
[`docs/IMPROVEMENTS.md`](docs/IMPROVEMENTS.md).

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
