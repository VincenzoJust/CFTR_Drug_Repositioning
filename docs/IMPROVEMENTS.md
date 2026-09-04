# What changed in the QSAR pipeline, and why

Short version: the training data had three quiet bugs that were deleting the
approved CFTR drugs and biasing labels upward, and the old pipeline was only
ever measuring itself with a random split, which is too easy for this kind of
chemistry data. Fixing the data and testing the model honestly gives two
numbers instead of one: **R² = 0.758** on a random split (looks great, but
optimistic) and **R² = 0.640** on a scaffold split (the real one — see #5).
Quote 0.640; it's the number that reflects how the model will actually
perform on new candidates.

Every number below was produced by running the pipeline, not estimated.

---

## The data fixes

**1. Lipinski's rule was deleting the approved drugs from training.**
Lipinski's rule of five predicts oral bioavailability — it has nothing to do
with whether a bioactivity measurement is valid. It was being used as a
cleaning filter anyway, and three of the four approved CFTR drugs fail it:
Ivacaftor (cLogP 5.08), Tezacaftor (MW 521), Elexacaftor (MW 598). A model
whose job is to score CFTR drugs was training without three of them. Fixed by
dropping Lipinski as a training filter — it's now just a descriptive column on
the candidates, which is what a drug-likeness heuristic is for. Effect: ChEMBL
513 → 803 compounds, merged dataset 1086 → **1591**, and all four approved
drugs are now in the training set.

**2. Censored measurements were treated as exact numbers.**
A value recorded as `>10000 nM` means "weaker than that," and `<10 nM` means
"stronger than that" — neither is a real number. Both were being parsed as if
they were. Fixed by keeping only exact (`=`/`~`) values in ChEMBL and
BindingDB, and by pulling ChEMBL's own `data_validity_comment` quality flag
(previously ignored).

**3. Replicate aggregation used the max, not the median.**
ChEMBL kept the *maximum* pIC50 per compound, which always picks the most
optimistic replicate. Papyrus and BindingDB already used the median. Switched
everything to median. Lumacaftor is the clean example: its five ChEMBL EC50
records are 2600, 2600, 2570, 2600 and 80 nM — max gives pEC50 7.1, median
gives 5.59, and 5.59 is the one that matches reality. A companion rule: any
compound whose replicates disagree by more than 1 log unit, or where two
databases disagree by more than 1 log unit, is dropped as unreliable.

**4. The target isn't really a "pIC50."** The ChEMBL pull is 87% functional
EC50 (potentiator-style assays), with smaller amounts of IC50, Kd and Ki mixed
in from 518 different assays. The column is still called `pIC50` for
continuity, but it's a pooled *pActivity*. Two indicator features
(`is_EC50`, `is_IC50`) are now fed to the model so it can learn to offset the
difference between assay types instead of treating them as identical.

---

## Honest evaluation

**5. Random cross-validation was inflating the R².** Random k-fold scatters
close chemical relatives across folds, so the model gets tested on molecules
nearly identical to ones it already memorized. This dataset has only 812
distinct chemical scaffolds across 1591 compounds — lots of near-duplicates.
Fix: evaluate under three schemes — random KFold, a scaffold-based split
(molecules sharing a chemical skeleton always land on the same side), and a
Butina similarity-cluster split. The seven repositioning candidates are all
structurally unrelated to anything in training (max similarity 0.26–0.33), so
the scaffold split is the one that matches how the model is actually used —
not just a stricter number, the *right* one.

A y-scrambling check (shuffle the labels, retrain) gives R² ≈ −0.10 on the
scaffold split — confirms the signal is real, not a fluke of the features.

**6. Model comparison, with a rule fixed before looking at results.** Random
Forest, XGBoost, tuned HistGradientBoosting, an SVR with a Tanimoto kernel,
and a neural network (MLP) were all compared under all three splits.
Selection rule declared in advance: whichever wins on the **scaffold split**
becomes production — a model that only wins on the easy random split is
memorizing analogues, not learning chemistry.

| Model | Random | **Scaffold** | Butina |
|---|---|---|---|
| Random Forest | 0.746 | 0.624 | 0.578 |
| XGBoost | 0.756 | 0.634 | 0.592 |
| HistGradientBoosting | 0.736 | 0.625 | 0.565 |
| SVR (Tanimoto) | 0.714 | 0.579 | 0.496 |
| MLP (neural net) | 0.400 | 0.047 | −0.129 |
| **Ensemble (RF+XGB+SVR)** | **0.758** | **0.640** | 0.581 |

The ensemble (simple average of RF, XGBoost, and the SVR) wins and is now the
production model. The neural network is a clean negative result, not a bug:
seven architectures were tried, the best reaches only R² = 0.047 on the
scaffold split, because ~1600 compounds is nowhere near enough data for a
neural net to beat gradient boosting on fingerprint features. A separate
feature-engineering sweep (radius-3 fingerprints, 2048 bits, MACCS keys, the
full 215-descriptor RDKit block) also went nowhere — every variant landed
within one fold's noise of the original 1031-feature set, and the biggest
feature sets were actually worse. Feature *selection* was deliberately not
attempted, since selecting on the same CV you then report is a leak.

---

## What the model can and can't tell you

**7. The QSAR term can't rank the 7 candidates against each other.** Their
predicted activities span only 0.374 units, while the per-compound
uncertainty (spread across the Random Forest's trees) is 0.22–0.56 — bigger
than the gap between candidates. In plain terms: the differences between the
candidates are noise, not signal, and the final ranking is driven almost
entirely by the docking score. The 60/40 docking/QSAR weighting was
deliberately **not** re-tuned after seeing this — that would be fitting the
score to the answer. Instead every prediction now ships with a 90% interval,
shown as error bars on the ranking figure.

**8. Control validation is now per-molecule, and split by mechanism.** The
old check compared group averages, which can pass even if one individual
positive control scores below the negatives — fixed to check each molecule
individually. It also now separates controls by mechanism: Ivacaftor is a
*potentiator*, while Lumacaftor/Tezacaftor/Elexacaftor are *correctors*. Since
the training data is ~87% potentiator-style assay, correctors are expected to
score lower — that's not a model failure, it's the pooled-endpoint limitation
from #4. Only the potentiator is required to clear the negative controls, and
it does with a clear margin (Ivacaftor 6.51 vs. best negative 5.85).

---

## Reproducibility and environment

- **A broken conda BLAS build was silently killing Python** on any numpy
  matrix multiply (even 200×200) — this broke the neural-network comparison
  with no traceback. Fixed by switching to OpenBLAS.
- **Papyrus was re-downloading itself every run** because of a version-string
  comparison bug. Fixed to try reading the cache first and only download on
  failure.
- **The notebook is now cache-first, not network-first.** Source data is
  cached under `data/cache/`, with a `REFRESH_DATA` flag to force a
  re-download, so results are pinned to a fixed snapshot instead of shifting
  between runs.
- **Papyrus was downloading 7.6 GB it never uses** — descriptor files
  (mordred, CDDD, mold2, etc.) that this notebook doesn't open. Of 9.71 GB on
  disk, only 1.30 GB was ever read. Fixed to skip those; the CFTR-relevant
  slice actually used is 0.25 MB and is committed to the repo, so a fresh
  clone runs the whole notebook offline with no large download.
- **Provenance is now recorded.** Each cached source carries a version/hash
  sidecar (ChEMBL release, Papyrus version, a SHA256 of the BindingDB export,
  which has no version string of its own), collected into
  `results/data_provenance.json`.
- ChEMBL API outages no longer crash the notebook at import time — it falls
  back to the cached data.
- Figures now save uniformly to `results/figures/` at 300 dpi.

---

## Limitations that are still there, on purpose

- The `pIC50 ≥ 5` filter is kept, so the model never sees an inactive
  compound and can't predict below ~5 — negative controls land at 5.47–5.85
  because that's the floor of the training range, not because the model
  recognizes them as inactive.
- The pooled EC50/IC50/Kd/Ki endpoint (#4) is not separated by mechanism —
  there isn't enough IC50 data to model correctors and potentiators
  separately.
- Censored compounds are dropped rather than modeled (proper censored
  regression isn't supported by scikit-learn's Random Forest and was judged
  out of scope).
- All seven candidates are chemical extrapolations (max similarity 0.26–0.33
  to training data) — their predictions sit outside the model's reliable
  applicability domain.
- The uncertainty estimate comes from Random Forest tree variance, which
  tends to be under-dispersed (too narrow), not from a proper Bayesian model.
- Docking energies are hard-coded from an external AutoDock Vina run, and the
  ~0.3 kcal/mol differences between candidates are well inside Vina's own
  error margin (~2–3 kcal/mol).

None of the following was done, because each would have quietly inflated the
numbers: tuning hyperparameters, features, or the 60/40 weights while looking
at the candidate predictions; reporting only the random-split R² once the
scaffold number was known; applying the replicate-disagreement cutoff after
seeing that it helped; any change aimed specifically at making Lumacaftor
predict higher; or selecting features on the same CV that's later quoted as
the performance estimate.
