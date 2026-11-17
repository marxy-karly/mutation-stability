
import json
import time
from pathlib import Path

import joblib
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import (
    ExtraTreesRegressor,
    RandomForestRegressor,
    HistGradientBoostingRegressor,
)
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.dummy import DummyRegressor


# --------------------------------------------------
# Configuration
# --------------------------------------------------

SPLIT_DIR = Path("data/splits")
RESULT_DIR = Path("results/comparison")
RESULT_DIR.mkdir(parents=True, exist_ok=True)

TARGET = "ddG_ML"

FEATURES = {
    "F1_baseline": [
        "aa_original",
        "aa_mutated",
        "position",
        "sequence_length",
    ],
    "F2_sequence_context": [
        "aa_original",
        "aa_mutated",
        "position",
        "sequence_length",
        "residue_minus_2",
        "residue_minus_1",
        "residue_plus_1",
        "residue_plus_2",
    ],
}

CATEGORICAL = [
    "aa_original",
    "aa_mutated",
    "residue_minus_2",
    "residue_minus_1",
    "residue_plus_1",
    "residue_plus_2",
]

NUMERICAL = ["position", "sequence_length"]


# --------------------------------------------------
# Load identical observation splits
# --------------------------------------------------

train = pd.read_csv(
    SPLIT_DIR / "train_features.csv",
    low_memory=False,
    dtype={"WT_cluster": "string"},
)

validation = pd.read_csv(
    SPLIT_DIR / "validation_features.csv",
    low_memory=False,
    dtype={"WT_cluster": "string"},
)

test = pd.read_csv(
    SPLIT_DIR / "test_features.csv",
    low_memory=False,
    dtype={"WT_cluster": "string"},
)

print("Dataset sizes:")
print("Train:", len(train))
print("Validation:", len(validation))
print("Test:", len(test))


# --------------------------------------------------
# Models
# --------------------------------------------------

def get_models():
    return {
        "Ridge": Ridge(alpha=10),

        "RandomForest": RandomForestRegressor(
            n_estimators=100,
            min_samples_leaf=2,
            n_jobs=-1,
            random_state=42,
        ),

        "ExtraTrees": ExtraTreesRegressor(
            n_estimators=100,
            min_samples_leaf=2,
            n_jobs=-1,
            random_state=42,
        ),

        "HistGradientBoosting": HistGradientBoostingRegressor(
            max_iter=100,
            max_leaf_nodes=31,
            l2_regularization=1.0,
            early_stopping=True,
            random_state=42,
        ),
    }


def build_pipeline(features, model):
    categorical = [f for f in features if f in CATEGORICAL]
    numerical = [f for f in features if f in NUMERICAL]

    preprocessing = ColumnTransformer([
        (
            "categorical",
            OneHotEncoder(
                handle_unknown="ignore",
                sparse_output=False,
            ),
            categorical,
        ),
        (
            "numerical",
            Pipeline([
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
            ]),
            numerical,
        ),
    ])

    return Pipeline([
        ("preprocessing", preprocessing),
        ("model", model),
    ])


def calculate_metrics(y_true, predictions):
    return {
        "MAE": mean_absolute_error(y_true, predictions),
        "RMSE": mean_squared_error(y_true, predictions) ** 0.5,
        "R2": r2_score(y_true, predictions),
    }


# --------------------------------------------------
# Train and evaluate on validation data
# --------------------------------------------------

results = []
best_score = float("inf")
best_pipeline = None
best_name = None
best_features = None

y_train = train[TARGET]
y_validation = validation[TARGET]

for feature_name, features in FEATURES.items():

    X_train = train[features]
    X_validation = validation[features]

    for model_name, model in get_models().items():

        print(f"\nTraining {model_name} | {feature_name}")
        start = time.time()

        pipeline = build_pipeline(features, model)
        pipeline.fit(X_train, y_train)

        predictions = pipeline.predict(X_validation)
        metrics = calculate_metrics(y_validation, predictions)

        elapsed = time.time() - start

        result = {
            "Feature_Set": feature_name,
            "Model": model_name,
            **metrics,
            "Training_Time_Seconds": round(elapsed, 2),
        }

        results.append(result)

        print("MAE:", round(metrics["MAE"], 5))
        print("RMSE:", round(metrics["RMSE"], 5))
        print("R2:", round(metrics["R2"], 5))
        print("Time:", round(elapsed, 2), "seconds")

        if metrics["MAE"] < best_score:
            best_score = metrics["MAE"]
            best_pipeline = pipeline
            best_name = model_name
            best_features = feature_name


# --------------------------------------------------
# Save validation results
# --------------------------------------------------

results_df = pd.DataFrame(results)
results_df = results_df.sort_values("MAE")

results_df.to_csv(
    RESULT_DIR / "validation_comparison.csv",
    index=False,
)

print("\nVALIDATION RANKING")
print(results_df.to_string(index=False))

print("\nBest model:", best_name)
print("Best feature set:", best_features)


# --------------------------------------------------
# Evaluate selected model on untouched test data
# --------------------------------------------------

selected_features = FEATURES[best_features]

X_test = test[selected_features]
y_test = test[TARGET]

test_predictions = best_pipeline.predict(X_test)
test_metrics = calculate_metrics(y_test, test_predictions)

test_results = {
    "Selected_Model": best_name,
    "Selected_Feature_Set": best_features,
    "Test_Observations": len(test),
    **test_metrics,
}

with open(RESULT_DIR / "test_results.json", "w") as f:
    json.dump(test_results, f, indent=4)

joblib.dump(
    best_pipeline,
    RESULT_DIR / "best_model.joblib",
)

print("\nFINAL TEST RESULTS")
print(json.dumps(test_results, indent=4))

print("\nSaved results to:", RESULT_DIR)

