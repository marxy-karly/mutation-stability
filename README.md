# Predicting the stability change of single mutations in unseen protein families

ES 335 Machine Learning course project.

**Question.** Given a protein's sequence and one amino-acid substitution in it, how well can the change in folding
stability (ΔΔG, `ddG_ML`, kcal/mol, positive = more stable) be predicted for protein families the model never saw?

- **Website:** `docs/` (7 pages: overview, guided tour, split, errors, compare, results, reproduce)
- **Live site:** _add your GitHub Pages URL_

## Repository layout

| Path | What it is |
|---|---|
| `src/preprocess.py` | Keeps single substitutions (`^[A-Z]\d+[A-Z]$`) with a numeric `ddG_ML`; extracts residues, position, length; validity checks |
| `src/split_data.py` | Two-stage `GroupShuffleSplit` on `WT_cluster`, seed 42 → train / validation / test |
| `src/train_models.py` | First baselines on F1 (Dummy, Ridge, Extra Trees) |
| `src/add_sequence_features.py` | Adds two residues on each side of the mutation (F2), `EDGE` past the chain ends |
| `src/update_splits.py` | Carries the same split assignments to the F2 files |
| `src/compare_models.py` | Four algorithms × F1/F2 on validation |
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
- The CSVs are not committed (GitHub rejects files over 100 MB). Download the zip and run the scripts below.

## Setup and run

```bash
pip install -r requirements.txt
python src/preprocess.py
python src/split_data.py
python src/train_models.py
python src/add_sequence_features.py
python src/update_splits.py
python src/compare_models.py
python src/evaluate.py               # writes results/ and docs/results.js
python -m http.server -d docs        # open http://localhost:8000
node tests/check_js.mjs              # optional
```

## What `src/evaluate.py` does

The split is never changed. Model selection uses validation MAE only; the test set is only scored.

1. **Baselines on validation and test:** mean; mean of the introduced residue; substitution table (mean ΔΔG of each original → new pair).
2. **Models:** Ridge, decision tree, random forest, extra trees, AdaBoost, HistGradientBoosting, and XGBoost if installed, each on F1 and F2. One fixed setting each (in `algorithms()`).
3. **Uncertainty:** 95% intervals from resampling whole test clusters (2,000 draws), for each score and for each model's difference from the mean and from the substitution table.
4. **Leakage check:** 5-fold CV on training + validation rows, grouped by cluster vs random.
5. **Repeated splits:** 5 grouped splits of training + validation (test untouched).
6. **Error analysis:** test errors by substitution type, designed vs natural, effect size, and cluster; permutation importance on validation.

## Results

See the website. Headline, scored on held-out clusters: _fill in after running `src/evaluate.py` on the full data_.

## Limitations

Simple sequence features with no 3D structure; only mutations with a measurable ΔΔG; few independent test clusters;
no hyperparameter tuning; the test set was scored once before for an earlier model; inferred column meanings.

## AI-use disclosure

Claude (Anthropic) wrote `src/evaluate.py`, the website (`docs/`, `tools/build_pages.py`) and `tests/check_js.mjs`,
and helped draft this README. I wrote the data pipeline in `src/preprocess.py` through `src/compare_models.py`,
chose the question and the split, and checked the results.
_Edit this to match exactly what you did._

## Citations

- Tsuboyama, K. et al. Mega-scale experimental analysis of protein folding stability in biology and design. *Nature* 620, 434–444 (2023). https://doi.org/10.1038/s41586-023-06328-6
- Data: https://doi.org/10.5281/zenodo.7401275 (CC BY 4.0)
- Pedregosa, F. et al. Scikit-learn: Machine Learning in Python. *JMLR* 12, 2825–2830 (2011).
