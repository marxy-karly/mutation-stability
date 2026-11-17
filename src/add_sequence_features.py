
import pandas as pd
from pathlib import Path

INPUT = Path("data/processed_dataset.csv")
OUTPUT = Path("data/processed_dataset_features.csv")

df = pd.read_csv(INPUT, low_memory=False)

def get_context(row):
    sequence = row["aa_seq"]
    position = int(row["position"]) - 1

    return pd.Series({
        "residue_minus_2": sequence[position - 2] if position >= 2 else "EDGE",
        "residue_minus_1": sequence[position - 1] if position >= 1 else "EDGE",
        "residue_plus_1": sequence[position + 1] if position + 1 < len(sequence) else "EDGE",
        "residue_plus_2": sequence[position + 2] if position + 2 < len(sequence) else "EDGE",
    })

context = df.apply(get_context, axis=1)
df = pd.concat([df, context], axis=1)

df.to_csv(OUTPUT, index=False)

print("Observations:", len(df))
print("New features:")
print(context.columns.tolist())
print("\nFeature value counts:")
for column in context.columns:
    print(f"\n{column}")
    print(context[column].value_counts().head())
print("\nSaved to:", OUTPUT)

