# Methodological improvements and limitations

A record of what was changed in the QSAR pipeline, **why**, and what each change
measurably did. Every number here was produced by running the pipeline, not
estimated. It also records what was tried and *didn't* work, which is the more
useful half.

Most of this maps directly onto a thesis "methodological improvements and
limitations" chapter.

---

## Headline result

| Metric | Before | After |
|---|---|---|
| Compounds in the training set | 1086 | **1591** |
| Approved CFTR drugs in training | 1 of 4 | **4 of 4** |
| R² — random 5-fold CV | 0.649 | **0.758** |
| R² — **scaffold split** (the honest one) | *not measured* | **0.640** |
| R² — Butina cluster split | *not measured* | 0.581 |
| Per-prediction uncertainty | *not reported* | reported everywhere |

The scaffold-split R² of 0.640 is the number to quote. It is *lower* than a
random-CV number would be, and that is the point: it is the one that estimates
error on genuinely new chemistry, which is what drug repositioning does.

---

## 1. Lipinski was deleting the approved drugs from the training set

**The problem.** `passes_lipinski()` was applied as a filter on the *training
data* in all three cleaning cells. Lipinski's rule of five predicts oral
bioavailability — it says nothing about whether a measurement is valid.

Measured on the ChEMBL pull (1002 unique compounds in nM):

| Filter set | Unique compounds |
|---|---|
| No filter | 1002 |
| `pIC50 >= 5` only | 897 |
| Lipinski only | 597 |
| Both (the old pipeline) | 513 |

Lipinski was **4× more destructive than the activity filter**. And critically:

| Drug | MW | cLogP | Passes Lipinski |
|---|---|---|---|
| Ivacaftor | 392 | **5.08** | No |
| Tezacaftor | **521** | 3.40 | No |
| Elexacaftor | **598** | 4.02 | No |
| Lumacaftor | 452 | 4.75 | Yes |

Three of the four approved CFTR drugs were being deleted from the training set
of a model whose purpose is to find CFTR modulators.

**The fix.** Removed from the three cleaning cells; retained as a descriptive
column on the candidates, which is what a drug-likeness heuristic is for.

**Effect.** ChEMBL 513 → 803 compounds; merged dataset 1086 → 1591. All four
approved drugs now present (the notebook asserts this explicitly).

## 2. Censored measurements were treated as exact

A value recorded as `>10000 nM` means the compound is *weaker* than that;
`<10 nM` means *stronger*. Both were being parsed as exact numbers — in ChEMBL
by ignoring `standard_relation` entirely, and in BindingDB by literally
stripping the `<`/`>` characters.

**The fix.** Fetch `standard_relation`, keep only `=` and `~`, and report the
cost. Also honour `data_validity_comment`, ChEMBL's own quality flag, which was
being ignored.

Tobit/censored regression was considered and rejected: it is not supported by
sklearn's RF and is out of scope. One honest sentence naming the excluded
compounds is worth more than a half-implemented survival model.

## 3. Aggregation was inconsistent and biased upward

ChEMBL kept the **maximum** pIC50 per compound (`sort_values(desc)` +
`drop_duplicates`), while Papyrus and BindingDB used the median. Max-aggregation
picks the most optimistic replicate.

Lumacaftor is the clean illustration: its five ChEMBL EC50 records are
2600, 2600, 2570, 2600 and **80** nM. Max gives pEC50 7.1; median gives 5.59.

**The fix.** Median everywhere, and Papyrus switched from `pchembl_value_Mean`
to `pchembl_value_Median`. Replicate spread (`pIC50_std`) and count (`n_meas`)
are now carried through to `dataset_full.csv`.

**A-priori quality rule** (threshold declared *before* looking at any model
result): drop compounds whose replicates disagree by more than 1 log unit, and
compounds where two databases disagree by more than 1 log unit.

## 4. The target is not really a pIC50

The ChEMBL pull is **1386 EC50, 158 IC50, 19 Kd, 3 Ki** across **518 distinct
assays** — about 87% functional EC50 from potentiator-style assays. The column
name `pIC50` is kept for continuity but the notebook now states plainly that it
is a pooled **pActivity**, and feeds `is_EC50` / `is_IC50` indicator features to
the model so it can offset the systematic difference between endpoints.

Separate corrector/potentiator models were considered and rejected: 158 IC50
compounds is not enough statistical power. The limitation is documented instead.

---

## 5. The R² was inflated by analogue leakage

**The problem.** Random *k*-fold CV scatters close analogues of the same chemical
series across folds, so the model is tested on molecules nearly identical to ones
it memorised. The 1591 compounds share only **812 Bemis-Murcko scaffolds** (593
singletons, largest class 36).

**The fix.** Every model is evaluated under three schemes — random KFold,
Bemis-Murcko scaffold `GroupKFold`, and Butina cluster `GroupKFold` (Tanimoto
0.65) — and the scaffold number is reported as the model's performance.

**Why scaffold CV is the *matching* regime here, not merely a stricter one:**
all seven candidates are structurally unrelated to the training data (max
Tanimoto 0.26–0.33). Scaffold CV measures exactly the error we will incur on
them. This turns a lower number into a methodological strength.

**y-scrambling check:** shuffle the labels, retrain, scaffold CV → R² = **−0.101**
(range −0.114 to −0.085). The signal is real, not an artefact of the features.

---

## 6. Model comparison — and a genuine negative result

All under 5-fold CV on 1591 compounds. **Selection rule declared before running:
best scaffold-split R² becomes production.**

| Model | Random KFold | **Scaffold** | Butina |
|---|---|---|---|
| Random Forest | 0.746 | 0.624 | 0.578 |
| XGBoost | 0.756 | 0.634 | 0.592 |
| HistGradientBoosting (tuned) | 0.736 | 0.625 | 0.565 |
| SVR (Tanimoto kernel) | 0.714 | 0.579 | 0.496 |
| **MLP (neural network)** | 0.400 | **0.047** | −0.129 |
| **Ensemble (RF+XGB+SVR)** | **0.758** | **0.640** | 0.581 |

Notes:

- The old notebook compared a **tuned** RF against a **default** HistGradientBoosting,
  which was not a fair test. Properly tuned, HGB moves 0.439 → 0.443 on the old
  data — it is simply the wrong family for 1024 sparse columns.
- XGBoost beats RF by ~0.01 — inside fold noise. Kept as a second family, not a
  headline.
- The **3-model ensemble** wins and becomes the production model.

### The neural network: a fair test, and a clear negative

Seven MLP configurations were tried under the scaffold split so the result could
not be dismissed as a strawman:

| Input / architecture | Scaffold R² |
|---|---|
| raw counts (512,128) | −0.200 |
| log1p (512,128) | −0.017 |
| **binary (512,128)** ← best | **+0.047** |
| log1p (256,) | −2.132 |
| log1p (128,64) | −0.361 |
| log1p, strong reg (256,) | −1.505 |
| binary, strong reg (256,) | −1.729 |

The best neural network reaches **0.047** against 0.624 for a plain Random
Forest. With ~1600 compounds this is expected — deep models generally need
10⁴–10⁶ examples to overtake gradient boosting on fingerprint input — and it is
reinforced by the feature ablation below: if richer representations buy nothing,
the ceiling is set by data quantity and label noise, not by architecture.

**This is a legitimate finding, not a failed experiment.** Transfer learning from
a pretrained chemical language model (e.g. frozen ChemBERTa embeddings + a light
head) is the one deep-learning route with real potential here, because the
chemistry knowledge would come from pretraining rather than from these 1591
compounds. It was scoped but not implemented.

---

## 7. What did NOT work: feature engineering

Tested under scaffold `GroupKFold`, RF(300, min_samples_leaf=2):

| Feature set | Dim | Scaffold R² |
|---|---|---|
| 7 desc + Morgan **bit** r2/1024 (original) | 1031 | 0.474 ± 0.052 |
| 7 desc + Morgan **count** r2/1024 ← adopted | 1031 | **0.478 ± 0.045** |
| 7 desc + count r3/2048 | 2055 | 0.466 ± 0.082 |
| full RDKit descriptors (215) + Morgan bit | 1239 | 0.461 ± 0.064 |
| full desc + count r3 + MACCS | 2430 | 0.452 ± 0.088 |
| full RDKit descriptors only | 215 | 0.414 ± 0.093 |

Every variant sits within one fold-SD of the baseline, and the **largest feature
sets are worse** with higher variance. The plan originally called feature
expansion the "biggest lever"; the measurement disproved that and it was dropped.

Only change adopted: binary → **count** fingerprints, worth it for the variance
reduction (0.052 → 0.045) rather than the +0.004.

Feature *selection* was deliberately not done: selecting on full-dataset CV and
then quoting that same CV is a leak.

---

## 8. The QSAR term cannot rank the candidates

The most consequential finding.

- Predicted activity across the 7 candidates spans **0.374** units (5.43–5.80).
- Per-compound uncertainty (RF tree SD) is **0.22–0.56**.

The candidates differ by **less than the uncertainty on any one of them**. The
40% QSAR weight in the integrated score is therefore ranking noise, and the
ordering is driven almost entirely by docking — which the pre-existing weight
sensitivity table independently confirms (the top of the ranking is unchanged at
`w_dock = 1.0`).

**Response.** The 60/40 weights were deliberately **not** re-tuned — doing so
after seeing the predictions would be fitting the score function to the answer.
Instead every prediction now carries a 90% interval, the ranking figure shows
error bars, and the notebook states in plain language that the candidate
differences are not significant.

## 9. Control validation was hiding a failure

The old check compared **group means**, which passes even when an individual
positive control scores below the negatives. It is now per molecule, against the
best negative control.

It also records which controls are **in the training set** — after fix #1 all
four approved drugs are, so their predictions are consistency checks rather than
independent validation. That is stated rather than glossed over.

### The Lumacaftor investigation

In the original results Lumacaftor (an approved corrector, Tanimoto 1.0) was
predicted 5.641 — *below* all three negative controls. This initially looked like
a model failure caused by the activity filter. **It is not.**

Lumacaftor's own ChEMBL records give a median **pEC50 of 5.59**; its label in
`dataset_full.csv` is 5.585 and the model predicted 5.64. The model was
reproducing its training label to within 0.06 log units — behaving perfectly.

The real issue is **endpoint pooling**. Lumacaftor is a *corrector* (it rescues
ΔF508 trafficking) being scored on an endpoint that is ~87% *potentiator*-style
functional EC50, where it is genuinely weak. Any change aimed at raising its
prediction would have been fitting to a control.

The notebook now splits positive controls by mechanism and requires only the
**potentiator** to clear the negatives. Current result: Ivacaftor 6.51 vs best
negative 5.85 — a margin of +0.66, and Tezacaftor (6.04) and Elexacaftor (6.87)
also clear it now that they are in the training set.

---

## 10. Environment and reproducibility fixes

- **Broken BLAS.** The conda env resolved to `mkl 2026.1.0`, under which *every*
  `numpy` matrix multiply — even 200×200 — silently killed the Python
  interpreter (exit 127, no traceback). This made the MLP, PCA and much of
  sklearn unusable. Fixed by switching the env to OpenBLAS
  (`conda install "libblas=*=*openblas"`).
- **Papyrus re-downloaded every run.** The guard was
  `if '05.7' not in get_downloaded_versions()`, but that helper returns
  `PapyrusVersion` objects whose string form is the release date (`2024.09.2`),
  so the comparison never matched and several GB were re-fetched each time. Now
  it attempts the read and downloads only on failure.
- **`papyrus-scripts` 3.0 API change.** `read_papyrus(chunksize=...)` no longer
  returns an iterable of chunks; any non-`None` value just switches to a single
  lazy polars `LazyFrame`. The old `for chunk in ...` loop raised
  `TypeError: LazyFrame is not subscriptable`.
- **ChEMBL outages no longer block the notebook.** Importing `new_client`
  fetches the API schema over the network, so an EBI HTTP 500 used to kill the
  notebook at the *imports* cell. The import now degrades gracefully and the
  download cell falls back to the cached `results/chembl_cftr_raw.csv`.
- **`Descriptors.CalcMolDescriptors()` hard-crashes** this RDKit build (the
  culprit is `AvgIpc`) — noted in the code so nobody reaches for it.
- Figures now go to `results/figures/` at a uniform 300 dpi via a `save_fig()`
  helper; previously they landed beside the CSVs at mixed 200/300 dpi.

---

## Open limitations

1. **The `pIC50 >= 5` filter is retained** (a deliberate decision to preserve the
   established methodology). The model therefore never sees an inactive compound
   and cannot predict below ~5. Negative controls land at 5.47–5.85 not because
   the model recognises them as inactive, but because that is the floor of its
   training range.
2. **Endpoint pooling.** EC50, IC50, Kd and Ki are pooled into one target. The
   indicator features help but do not separate correctors from potentiators;
   there is not enough IC50 data to model them separately.
3. **Censored compounds are discarded** rather than modelled.
4. **All candidates are extrapolations** — max Tanimoto 0.26–0.33 to anything in
   training. Their predictions are outside the reliable applicability domain.
5. **Uncertainty is approximate.** It comes from Random Forest tree variance,
   which is generally under-dispersed, and is centred on the ensemble prediction.
6. **Docking energies are hard-coded** from an external AutoDock Vina run, and
   ~0.3 kcal/mol differences are well inside Vina's error (~2–3 kcal/mol).

## Things deliberately not done

Because they would inflate the numbers dishonestly:

- Tuning hyperparameters, features or the 60/40 weights while looking at the
  control predictions or the candidate ranking.
- Reporting only the random-KFold R² once the scaffold number was known.
- Applying the `pIC50_std > 1` cut *after* seeing that it helped.
- Any change aimed at making Lumacaftor predict high.
- Selecting features on full-dataset CV and quoting that same CV as the estimate.
