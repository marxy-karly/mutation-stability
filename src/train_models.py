
import pandas as pd
import numpy as np
import json
from pathlib import Path

from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import Ridge
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    r2_score
)

# Configuration
DATA = Path("data/splits")
RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)

FEATURES = [
    "aa_original",
    "aa_mutated",
    "position",
    "sequence_length"
]

TARGET = "ddG_ML"

# Load partitions
train = pd.read_csv(DATA / "train.csv")
val = pd.read_csv(DATA / "validation.csv")
test = pd.read_csv(DATA / "test.csv")

X_train, y_train = train[FEATURES], train[TARGET]
X_val, y_val = val[FEATURES], val[TARGET]
X_test, y_test = test[FEATURES], test[TARGET]

categorical = ["aa_original", "aa_mutated"]
numerical = ["position", "sequence_length"]

# Preprocessing
preprocessor = ColumnTransformer([
    ("cat", OneHotEncoder(handle_unknown="ignore"), categorical),
    ("num", StandardScaler(), numerical)
])

# Models
models = {
    "Dummy": DummyRegressor(strategy="mean"),

    "Ridge": Pipeline([
        ("preprocess", preprocessor),
        ("model", Ridge(alpha=10))
    ]),

    "ExtraTrees": Pipeline([
        ("preprocess", preprocessor),
        ("model", ExtraTreesRegressor(
            n_estimators=100,
            min_samples_leaf=2,
            n_jobs=-1,
            random_state=42
        ))
    ])
}

def evaluate(model, X, y):
    pred = model.predict(X)

    return {
        "MAE": mean_absolute_error(y, pred),
        "RMSE": np.sqrt(mean_squared_error(y, pred)),
        "R2": r2_score(y, pred)
    }

results = {}

# Train and validate
for name, model in models.items():
    print(f"\nTraining {name}...")

    model.fit(X_train, y_train)

    metrics = evaluate(model, X_val, y_val)
    results[name] = metrics

    print("Validation:", metrics)

# Compare models
comparison = pd.DataFrame(results).T
comparison.to_csv(RESULTS / "validation_metrics.csv")

print("\nVALIDATION COMPARISON")
print(comparison.to_string())

# Select best model by validation MAE
best_name = comparison["MAE"].idxmin()
best_model = models[best_name]

print("\nBest model:", best_name)

# Final evaluation on untouched test set
test_metrics = evaluate(best_model, X_test, y_test)

print("\nFINAL TEST RESULTS")
print(test_metrics)

# Save test metrics
with open(RESULTS / "test_metrics.json", "w") as f:
    json.dump(test_metrics, f, indent=4)

# Save selected model
import joblib

joblib.dump(
    best_model,
    RESULTS / "best_model.joblib"
)

print("\nModel and results saved.")

