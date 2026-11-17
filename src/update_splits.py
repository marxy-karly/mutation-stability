
import pandas as pd
from pathlib import Path

DATA = Path("data/processed_dataset_features.csv")
SPLIT_DIR = Path("data/splits")

df = pd.read_csv(DATA, low_memory=False, dtype={"WT_cluster": "string"})

# Existing split assignments
for split in ["train", "validation", "test"]:
    old = pd.read_csv(
        SPLIT_DIR / f"{split}.csv",
        low_memory=False,
        dtype={"WT_cluster": "string"}
    )

    # Match observations using their unique protein and mutation identifiers
    keys = ["WT_name", "mut_type"]
    merged = old[keys].merge(
        df,
        on=keys,
        how="left",
        validate="one_to_one",
        indicator=True
    )

    if not (merged["_merge"] == "both").all():
        raise ValueError(f"Unmatched observations in {split}")

    merged = merged.drop(columns="_merge")
    merged.to_csv(SPLIT_DIR / f"{split}_features.csv", index=False)

    print(f"{split}: {len(merged)} observations")

print("\nExpanded splits created successfully.")

