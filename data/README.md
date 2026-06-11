# Input data

Most of the data is downloaded automatically by the notebook. Only **BindingDB**
needs a manual download.

## 1. BindingDB (manual — required for the full dataset)

BindingDB has no API, so the file is exported by hand:

1. Go to <https://www.bindingdb.org>.
2. Search by target: **CFTR** (UniProt **P13569**), or use the advanced search.
3. Export the results as a **TSV** (tab-separated) file.
4. Save it in this folder as:

   ```
   data/bindingdb_cftr.tsv
   ```

The notebook reads these columns (others are ignored):
`Ligand SMILES`, `IC50 (nM)`, `Ki (nM)`, `Kd (nM)`, `EC50 (nM)`.

If the file is not present, the notebook prints a warning and continues with
**ChEMBL + Papyrus only** — it will not crash.

## 2. ChEMBL (automatic)

Downloaded via the ChEMBL web client, target **CHEMBL4051** (CFTR), assay types
IC50 / Ki / EC50, organism *Homo sapiens*. A raw copy is cached to
`results/chembl_cftr_raw.csv`.

## 3. Papyrus (automatic)

Version **05.7** is downloaded and cached by `papyrus-scripts` on first run
(several GB). Filtered by UniProt **P13569**, high-quality records only.

## 4. Docking / structural inputs (external, not required to run the notebook)

The binding energies in section 11 come from **AutoDock Vina** runs performed
outside this notebook on:

- **PDB 6MSM** — apo wild-type CFTR (`dock_wt`)
- an **AlphaFold2 ΔF508** model (`dock_delta`)

These values, and the **PLIP** interaction counts in section 16, are hard-coded
in the notebook, so you do **not** need the structure files to reproduce the
final ranking.
