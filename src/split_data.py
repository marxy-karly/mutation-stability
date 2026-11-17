
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from pathlib import Path

df = pd.read_csv(
    "data/processed_dataset.csv",
    dtype={"WT_cluster": "string"},
    low_memory=False
)

# Remove whitespace from cluster identifiers
df["WT_cluster"] = df["WT_cluster"].str.strip()

# Ensure groups are valid
assert df["WT_cluster"].notna().all()

groups = df["WT_cluster"]

# 70% train+validation, 30% test
split1 = GroupShuffleSplit(
    n_splits=1,
    test_size=0.30,
    random_state=42
)

train_val_idx, test_idx = next(
    split1.split(df, groups=groups)
)

# 17.65% of remaining = approximately 15% overall
split2 = GroupShuffleSplit(
    n_splits=1,
    test_size=0.1765,
    random_state=42
)

train_rel, val_rel = next(
    split2.split(
        df.iloc[train_val_idx],
        groups=groups.iloc[train_val_idx]
    )
)

train_idx = train_val_idx[train_rel]
val_idx = train_val_idx[val_rel]

# Save datasets
Path("data/splits").mkdir(
    parents=True, exist_ok=True
)

for name, idx in [
    ("train", train_idx),
    ("validation", val_idx),
    ("test", test_idx)
]:
    part = df.iloc[idx]
    part.to_csv(
        f"data/splits/{name}.csv",
        index=False
    )

    print(
        f"{name}: {len(part):,} observations, "
        f"{part['WT_cluster'].nunique()} clusters"
    )

# Verify cluster isolation
sets = [
    set(df.iloc[idx]["WT_cluster"])
    for idx in [train_idx, val_idx, test_idx]
]

assert not sets[0] & sets[1]
assert not sets[0] & sets[2]
assert not sets[1] & sets[2]

print("\nCluster leakage check: PASSED")

