# Predicting the stability change of single mutations in unseen protein families

ES 335 Machine Learning course project.

**Question.** Given a protein's sequence and one amino-acid substitution in it, how well can the change in folding
stability (ΔΔG, `ddG_ML`, kcal/mol, positive = more stable) be predicted for protein families the model never saw?

- **Website:** `docs/` (7 pages: overview, guided tour, split, errors, compare, results, reproduce)
- **Live site:** https://marxy-karly.github.io/mutation-stability/

## Repository layout

| Path | What it is |
|---|---|
| `src/preprocess.py` | Keeps single substitutions (`^[A-Z]\d+[A-Z]$`) with a numeric `ddG_ML`; extracts residues, position, length; validity checks |
| `src/split_data.py` | Two-stage `GroupShuffleSplit` on `WT_cluster`, seed 42 → train / validation / test |
| `src/train_models.py` | First baselines on F1 (Dummy, Ridge, Extra Trees) |
| `src/add_sequence_features.py` | Adds two residues on each side of the mutation (F2), `EDGE` past the chain ends |
| `src/update_splits.py` | Carries the same split assignments to the F2 files |
| `src/compare_models.py` | Four algorithms × F1/F2 on validation |
| `src/features.py` | F3 (sequence-derived physicochemical and neighbourhood features) and F4 (burial in the AlphaFold model); `python src/features.py` reports structure coverage |
| `src/evaluate.py` | Final evaluation: everything shown on the website (see below) |
| `config.json` | Column names, file paths, seeds, bootstrap and repeat counts, author fields |
| `results/evaluation.json` | Every number on the website |
| `results/test_predictions.csv.gz` | Every test mutation with every model's prediction |
| `docs/` | The website. `docs/results.js` is written by `src/evaluate.py` |
| `tools/build_pages.py` | Regenerates the HTML pages after text edits |
| `tests/check_js.mjs` | Checks that the website's in-browser ridge matches scikit-learn |

## Data

- Tsuboyama et al. (2023), Zenodo record 7401275: https://zenodo.org/records/7401275
- File: `K50_dG_Dataset1_Dataset2.csv` inside `Processed_K50_dG_datasets.zip` (CC BY 4.0; the record asks users to register their use).
- 851,552 rows → 542,542 single substitutions → 376,918 with a numeric `ddG_ML` (165,624 hold `-`).
- Structures for F4: `AlphaFold_model_PDBs.zip` from the same record (13.7 MB), unzipped into `data/structures/`. They are matched to proteins by sequence, so file names do not matter. Without them F4 is skipped and everything else runs.
- The CSVs are not committed (GitHub rejects files over 100 MB). Download the zips and run the scripts below.

## Setup and run

```bash
pip install -r requirements.txt
python src/preprocess.py
python src/split_data.py
python src/train_models.py
python src/add_sequence_features.py
python src/update_splits.py
python src/compare_models.py
unzip AlphaFold_model_PDBs.zip -d data/structures   # structures for F4
python src/features.py               # optional: how many proteins have a structure
python src/evaluate.py               # writes results/ and docs/results.js
python -m http.server -d docs        # open http://localhost:8000
node tests/check_js.mjs              # optional
```

## What `src/evaluate.py` does

The split is never changed. Model selection uses validation MAE only; the test set is only scored.

1. **Baselines on validation and test:** mean; mean of the introduced residue; substitution table (mean ΔΔG of each original → new pair).
2. **Models:** Ridge, decision tree, random forest, extra trees, AdaBoost, HistGradientBoosting, and XGBoost if installed, each on F1 and F2. The algorithms in `extended_algorithms` (config.json) also run on F3 and F4. One fixed setting each (in `algorithms()`).
   - **F3** = F2 + how the swap changes hydrophobicity, volume, charge and helix propensity; relative position and distance to the chain end; hydrophobicity and charge of the neighbours, including the residues on the same face of a helix (i±3, i±4) or strand (i±2).
   - **F4** = F3 + at the mutated position in the AlphaFold model of the unmutated protein: residues within 8 and 10 Å, half-sphere exposure, relative solvent-accessible surface area (Shrake–Rupley), mean hydrophobicity of the 3D neighbours, distance from the protein centre, pLDDT, a helix/strand geometry measure, burial rank within the protein, and burial × (hydrophobicity, volume, helix-propensity) change. Missing values (no structure) are filled with the training median plus a 0/1 flag, inside the pipeline.
3. **Uncertainty:** 95% intervals from resampling whole test clusters (2,000 draws), for each score and for each model's difference from the mean and from the substitution table.
4. **Leakage check:** 5-fold CV on training + validation rows, grouped by cluster vs random.
5. **Repeated splits:** 5 grouped splits of training + validation (test untouched), including the selected algorithm on every feature set.
7. **Tuning and ensemble:** a small grid for HistGradientBoosting and XGBoost on the richest feature set, chosen on validation MAE (no early stopping, which would use a random rather than grouped slice); an ensemble averaging the tuned HistGradientBoosting, tuned XGBoost and random forest. Both compete in the same validation-based selection.
6. **Error analysis:** test errors by substitution type, designed vs natural, effect size, burial, and cluster; permutation importance on validation, per column and per kind of information.

## Results

Final run (`python src/evaluate.py` on the full data). Test MAE in kcal/mol, 95% intervals from resampling whole test clusters.

| Method | Validation MAE | Test MAE | 95% interval |
|---|---|---|---|
| Mean (baseline) | 0.714 | 0.847 | 0.745 – 0.915 |
| Substitution table (baseline) | 0.562 | 0.701 | 0.603 – 0.767 |
| XGBoost (F2), best before structure | 0.547 | 0.685 | 0.601 – 0.741 |
| HistGradientBoosting (F4) | 0.466 | 0.594 | 0.521 – 0.640 |
| **Ensemble (F4), selected by validation MAE** | **0.455** | **0.575** | **0.506 – 0.618** |

- The selected ensemble beats the substitution table by 0.126 kcal/mol (95% interval 0.084 to 0.157). Before the
  structure features (F4), no model reliably beat it.
- Burial in the AlphaFold models (F4) gave the largest gain; the sequence neighbours (F2) and sequence-derived numbers (F3)
  changed MAE by about 0.01 each; tuning and the ensemble added about 0.02.
- Leakage check: random folds make the ensemble look 22% better than folds that keep clusters together (0.415 against 0.529).
- 387 of 396 wild-type proteins matched an AlphaFold model; 96–100% of rows in each split have structure features.

Every number is in `results/evaluation.json`; every test prediction in `results/test_predictions.csv.gz`.

## Limitations and sources of prediction error

The final ensemble has a test MAE of 0.575 kcal/mol. The website (Results page) gives each point with the supporting numbers.

1. **Incomplete structural environment.** Burial, solvent accessibility and neighbourhood features describe how crowded and
   exposed a residue is, not which interactions it makes (hydrogen bonds, salt bridges, packing, disulfides, specific contacts).
2. **Wild-type, not mutant, structure.** Mutation-induced rearrangements (clashes, lost hydrogen bonds, changed exposure) are not visible.
3. **Uncertainty in the target.** `ddG_ML` is inferred from proteolysis measurements through a model of the assay; its noise
   floor was not estimated.
4. **Generalization to unseen clusters.** Validation MAE 0.455 but test MAE 0.575: the test families are less like the training data.
5. **Underprediction of extreme changes.** Errors grow with the size of the real effect; predictions are pulled toward the centre.
6. **Directions:** ESM-2 embeddings, wild-type vs mutant structure features, residue-interaction features, uncertainty-aware
   training, refitting the final model on training + validation.

The contribution of each to the remaining error has not been quantified separately. Other limitations: only mutations with a
measurable ΔΔG; few independent test clusters; uneven split proportions; small tuning grid for two algorithms; the test set was
scored after each round of feature work (choices were made on validation, and the repeated grouped splits never use the test set);
inferred column meanings.

## Author, contribution and AI-use text

The website reads `author`, `contribution`, `ai_use` and `repo_url` from `config.json`. After editing them, run
`python tools/build_pages.py` (a few seconds; no need to rerun the evaluation). Leave a blank line between paragraphs.

## AI-use disclosure

Generative AI was used for the following:
- Create guides to understand the research paper methods
- Acting as a methodological mentor to understand the underlying biological/physical principles.
- Build the website
- Assistance in writing `README.md`. Assistance in writing the documentation for this project.
- Assistance in feature engineering and understanding the biological relevance of features.
- Assistance in writing `src/` files.
## Citations

- Tsuboyama, K. et al. Mega-scale experimental analysis of protein folding stability in biology and design. *Nature* 620, 434–444 (2023). https://doi.org/10.1038/s41586-023-06328-6
- Data: https://doi.org/10.5281/zenodo.7401275 (CC BY 4.0)
- Jumper, J. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589 (2021).
- Kyte, J. & Doolittle, R. F. A simple method for displaying the hydropathic character of a protein. *J. Mol. Biol.* 157, 105–132 (1982).
- Zamyatnin, A. A. Protein volume in solution. *Prog. Biophys. Mol. Biol.* 24, 107–123 (1972).
- Pace, C. N. & Scholtz, J. M. A helix propensity scale based on experimental studies of peptides and proteins. *Biophys. J.* 75, 422–427 (1998).
- Hamelryck, T. An amino acid has two sides: a new 2D measure provides a different view of solvent exposure. *Proteins* 59, 38–48 (2005).
- Shrake, A. & Rupley, J. A. Environment and exposure to solvent of protein atoms. *J. Mol. Biol.* 79, 351–371 (1973).
- Tien, M. Z. et al. Maximum allowed solvent accessibilities of residues in proteins. *PLoS ONE* 8, e80635 (2013).
- Pedregosa, F. et al. Scikit-learn: Machine Learning in Python. *JMLR* 12, 2825–2830 (2011).
