
import pandas as pd
import numpy as np
from pathlib import Path

# Configuration
INPUT = Path(
    "~/Dataset/Processed_K50_dG_datasets/"
    "single_point_mutations.csv"
).expanduser()

OUTPUT = Path("data/processed_dataset.csv")

TARGET = "ddG_ML"

# Load only relevant columns
columns = [
    "mut_type",
    "aa_seq",
    "aa_seq_full",
    "WT_name",
    "WT_cluster",
    TARGET
]

df = pd.read_csv(
    INPUT,
    usecols=columns,
    low_memory=False,
    dtype={"WT_cluster": "string"}
)

print("Initial observations:", len(df))

# 1. Convert target to numeric
df[TARGET] = pd.to_numeric(
    df[TARGET], errors="coerce"
)

df = df.dropna(subset=[TARGET]).copy()

print("Valid targets:", len(df))

# 2. Parse mutation notation
mutation = df["mut_type"].str.extract(
    r"^([A-Z])(\d+)([A-Z])$"
)

df["aa_original"] = mutation[0]
df["position"] = pd.to_numeric(
    mutation[1], errors="coerce"
)
df["aa_mutated"] = mutation[2]

# 3. Validate mutation parsing
df = df.dropna(
    subset=["aa_original", "position", "aa_mutated"]
).copy()

# 4. Basic sequence features
df["sequence_length"] = (
    df["aa_seq"].str.len()
)

# 5. Remove rows missing essential information
df = df.dropna(
    subset=["WT_name", "WT_cluster", "aa_seq"]
).copy()

# 6. Check sequence validity
valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

df["valid_sequence"] = df["aa_seq"].apply(
    lambda s: (
        isinstance(s, str)
        and len(s) > 0
        and set(s.upper()).issubset(valid_amino_acids)
    )
)

invalid_count = (~df["valid_sequence"]).sum()

print("Invalid sequences:", invalid_count)

df = df[df["valid_sequence"]].copy()

# 7. Check mutation position bounds
df["position_valid"] = (
    (df["position"] >= 1)
    & (df["position"] <= df["sequence_length"])
)

print(
    "Positions outside sequence bounds:",
    (~df["position_valid"]).sum()
)

# Retain only observations with valid position bounds
df = df[df["position_valid"]].copy()

# 8. Save modelling dataset
df["WT_cluster"] = df["WT_cluster"].str.strip()

df = df.drop(
    columns=["valid_sequence", "position_valid"]
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(OUTPUT, index=False)

# 9. Report
print("\n========== PREPROCESSING REPORT ==========")
print("Final observations:", len(df))
print("Unique proteins:", df["WT_name"].nunique())
print("Unique clusters:", df["WT_cluster"].nunique())
print("Unique mutations:", df["mut_type"].nunique())
print("Missing targets:", df[TARGET].isna().sum())
print("Target mean:", df[TARGET].mean())
print("Target std:", df[TARGET].std())
print("Output:", OUTPUT)

