"""Final evaluation. Run after preprocess.py, split_data.py, add_sequence_features.py and update_splits.py.

    python src/evaluate.py

It reads the frozen splits in data/splits/*_features.csv (never re-splits the test set) and writes:
    results/evaluation.json          every number the website shows
    results/test_predictions.csv.gz  every test mutation with every model's prediction
    docs/results.js                  the same results, for the website

What it adds to compare_models.py:
    1. Baselines scored on the TEST set too: the mean, and a substitution table
       (average ddG of each original -> new amino-acid pair in the training rows).
    2. Decision tree, AdaBoost and (if installed) XGBoost next to the four algorithms already compared.
    3. A 95% interval for every test score, by resampling whole test clusters (cluster bootstrap).
    4. A leakage check: the same models under random folds and under folds that keep clusters together.
    5. Repeated grouped splits of the training + validation rows (the test set is not touched).
    6. Error analysis on the test set and permutation importance on the validation set.

Model selection is unchanged: the learned model with the lowest VALIDATION MAE. The test set is only scored.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import (AdaBoostRegressor, ExtraTreesRegressor, HistGradientBoostingRegressor,
                              RandomForestRegressor)
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, GroupShuffleSplit, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.tree import DecisionTreeRegressor

ROOT = Path(__file__).resolve().parent.parent
AAS = "ACDEFGHIKLMNPQRSTVWY"
CATS = list(AAS) + ["EDGE"]
CLASS = {**{a: "hydrophobic" for a in "AILMFVWC"}, **{a: "polar" for a in "STNQYH"}, "K": "positive", "R": "positive",
         "D": "negative", "E": "negative", "G": "special (G, P)", "P": "special (G, P)"}


# ----------------------------------------------------------------------------- models
def algorithms(seed):
    """The learned models. If compare_models.py used other settings, copy them here so the numbers match."""
    algs = {
        "Ridge": lambda: Ridge(alpha=1.0),
        "Decision tree": lambda: DecisionTreeRegressor(max_depth=10, min_samples_leaf=50, random_state=seed),
        "Random forest": lambda: RandomForestRegressor(n_estimators=100, min_samples_leaf=5, n_jobs=-1, random_state=seed),
        "Extra trees": lambda: ExtraTreesRegressor(n_estimators=100, min_samples_leaf=5, n_jobs=-1, random_state=seed),
        "AdaBoost": lambda: AdaBoostRegressor(n_estimators=100, learning_rate=0.1, random_state=seed),
        "HistGradientBoosting": lambda: HistGradientBoostingRegressor(random_state=seed),
    }
    try:
        from xgboost import XGBRegressor
        algs["XGBoost"] = lambda: XGBRegressor(n_estimators=300, max_depth=6, learning_rate=0.1, subsample=0.8,
                                               n_jobs=-1, random_state=seed)
    except ImportError:
        print("xgboost is not installed: skipping XGBoost (pip install xgboost to include it)")
    return algs


def pipeline(fs, model):
    """Same preprocessing as train_models.py: one-hot for amino acids, standard scaling for numbers.
    Both are fitted on the training rows only, because they sit inside the pipeline."""
    pre = ColumnTransformer([("cat", OneHotEncoder(handle_unknown="ignore", sparse_output=False), fs["cat"]),
                             ("num", StandardScaler(), fs["num"])])
    return Pipeline([("pre", pre), ("model", model)])


class GroupMean:
    """Baseline: the training average of y for each value of a key; unseen keys fall back to another baseline."""
    def __init__(self, key_fn, fallback=None):
        self.key_fn, self.fallback = key_fn, fallback

    def fit(self, df, y):
        self.mu = float(np.mean(y))
        self.table = pd.Series(np.asarray(y)).groupby(self.key_fn(df).values).mean()
        if self.fallback:
            self.fallback.fit(df, y)
        return self

    def predict(self, df):
        p = self.key_fn(df).map(self.table).to_numpy(dtype=float)
        fb = self.fallback.predict(df) if self.fallback else np.full(len(df), self.mu)
        return np.where(np.isnan(p), fb, p)


def baselines(c):
    mean = lambda: GroupMean(lambda d: pd.Series("all", index=d.index))
    mutant = lambda: GroupMean(lambda d: d[c["mutated"]], mean())
    table = lambda: GroupMean(lambda d: d[c["original"]] + ">" + d[c["mutated"]], mutant())
    return {"Mean": mean, "Mutant residue mean": mutant, "Substitution table": table}


# ----------------------------------------------------------------------------- data
def load(cfg):
    c = cfg["columns"]
    need = [c[k] for k in ("target", "group", "protein", "mutation", "original", "mutated", "position", "length")] + c["context"]
    parts = {}
    for part in ("train", "validation", "test"):
        p = ROOT / cfg["splits"][part]
        if not p.exists():
            sys.exit(f"Missing {p}. Run the earlier scripts in src/ first (see README).")
        df = pd.read_csv(p, low_memory=False)
        miss = [x for x in need if x not in df.columns]
        if miss:
            sys.exit(f"{p.name} has no column(s) {miss}.\nIts columns are: {list(df.columns)}\nEdit 'columns' in config.json.")
        for col in [c["group"], c["protein"], c["mutation"], c["original"], c["mutated"]] + c["context"]:
            df[col] = df[col].astype(str).str.strip()
        df[c["target"]] = pd.to_numeric(df[c["target"]], errors="raise")
        parts[part] = df.reset_index(drop=True)
    g = {k: set(v[c["group"]]) for k, v in parts.items()}
    shared = (g["train"] & g["validation"]) | (g["train"] & g["test"]) | (g["validation"] & g["test"])
    if shared:
        sys.exit(f"Leakage: clusters {sorted(shared)[:5]} appear in more than one split.")
    return parts


def mae(y, p):
    return float(np.mean(np.abs(np.asarray(y) - np.asarray(p))))


def rmse(y, p):
    return float(np.sqrt(np.mean((np.asarray(y) - np.asarray(p)) ** 2)))


def r2(y, p):
    y = np.asarray(y)
    return float(1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2))


def spear(y, p):
    return None if np.std(p) < 1e-12 else float(spearmanr(y, p)[0])


def within_protein_spearman(prot, y, p, min_rows=20):
    vals = []
    for _, i in pd.Series(range(len(y))).groupby(prot).groups.items():
        i = np.asarray(list(i))
        if len(i) >= min_rows and np.std(p[i]) > 1e-12 and np.std(y[i]) > 1e-12:
            vals.append(spearmanr(y[i], p[i])[0])
    return float(np.mean(vals)) if vals else None


def design_f1(df, c):
    """F1 as a plain matrix, columns in a fixed order. Used only to check the website's in-browser ridge."""
    cols = []
    for col in (c["original"], c["mutated"]):
        cols += [(df[col] == a).astype(float).to_numpy() for a in CATS]
    cols += [df[c["position"]].astype(float).to_numpy(), df[c["length"]].astype(float).to_numpy()]
    return np.column_stack(cols)


# ----------------------------------------------------------------------------- main
def main():
    t0 = time.time()
    cfg = json.loads((ROOT / "config.json").read_text(encoding="utf-8"))
    c, seed = cfg["columns"], cfg["seed"]
    D = load(cfg)
    tr, va, te = D["train"], D["validation"], D["test"]
    T = c["target"]
    ytr, yva, yte = tr[T].to_numpy(float), va[T].to_numpy(float), te[T].to_numpy(float)
    FS = {"F1": {"cat": [c["original"], c["mutated"]], "num": [c["position"], c["length"]]}}
    FS["F2"] = {"cat": FS["F1"]["cat"] + c["context"], "num": FS["F1"]["num"]}
    print(f"train {len(tr)}  validation {len(va)}  test {len(te)}")

    # ---- 1. fit every model on the training rows; score validation and test
    makers = {}
    for name, mk in baselines(c).items():
        makers[name] = ("baseline", "–", mk)
    for alg, mk in algorithms(seed).items():
        for fs in ("F1", "F2"):
            makers[f"{alg} ({fs})"] = ("model", fs, (lambda mk=mk, fs=fs: pipeline(FS[fs], mk())))
    pva, pte = {}, {}
    for name, (kind, fs, mk) in makers.items():
        t1 = time.time()
        m = mk().fit(tr, ytr)
        pva[name], pte[name] = np.asarray(m.predict(va), float), np.asarray(m.predict(te), float)
        print(f"  {name:32s} val MAE {mae(yva, pva[name]):.4f}  ({time.time() - t1:.0f}s)", flush=True)
    learned = [n for n in makers if makers[n][0] == "model"]
    selected = min(learned, key=lambda n: mae(yva, pva[n]))
    sel_alg, sel_fs = selected.rsplit(" (", 1)[0], makers[selected][1]
    REF = "Substitution table"
    print(f"selected by validation MAE: {selected}")

    # ---- 2. cluster bootstrap on the test set: resample whole clusters
    ug, inv = np.unique(te[c["group"]], return_inverse=True)
    cnt = np.bincount(inv).astype(float)
    W = np.random.default_rng(seed).multinomial(len(ug), np.full(len(ug), 1 / len(ug)), size=cfg["n_bootstrap"]).astype(float)
    den = W @ cnt
    boot = {n: W @ np.bincount(inv, weights=np.abs(yte - pte[n])) / den for n in makers}
    ci = lambda a: [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]
    prot = te[c["protein"]].to_numpy()
    models = []
    for n, (kind, fs, _) in makers.items():
        models.append({
            "name": n, "kind": kind, "features": fs, "algorithm": n.rsplit(" (", 1)[0],
            "val_mae": mae(yva, pva[n]), "val_rmse": rmse(yva, pva[n]), "val_r2": r2(yva, pva[n]),
            "test_mae": mae(yte, pte[n]), "test_mae_ci": ci(boot[n]), "test_rmse": rmse(yte, pte[n]), "test_r2": r2(yte, pte[n]),
            "test_spearman": spear(yte, pte[n]), "test_within_spearman": within_protein_spearman(prot, yte, pte[n]),
            "diff_vs_mean": mae(yte, pte[n]) - mae(yte, pte["Mean"]), "diff_vs_mean_ci": ci(boot[n] - boot["Mean"]),
            "diff_vs_ref": mae(yte, pte[n]) - mae(yte, pte[REF]), "diff_vs_ref_ci": ci(boot[n] - boot[REF]),
        })

    # ---- 3. leakage check on training + validation rows (test untouched)
    pool = pd.concat([tr, va], ignore_index=True)
    ypool = pool[T].to_numpy(float)
    rng = np.random.default_rng(seed)
    li = np.sort(rng.choice(len(pool), min(cfg["leak_rows"], len(pool)), replace=False))
    L, yL, gL = pool.iloc[li].reset_index(drop=True), ypool[li], pool[c["group"]].to_numpy()[li]
    k = cfg["leak_folds"]
    folds = {"grouped": list(GroupKFold(n_splits=k).split(L, groups=gL)),
             "random": list(KFold(n_splits=k, shuffle=True, random_state=seed).split(L))}
    leak_names = list(dict.fromkeys(["Mean", REF, "Ridge (F1)", selected]))
    leak = []
    for n in leak_names:
        row = {"name": n}
        for how, fl in folds.items():
            row[how] = float(np.mean([mae(yL[b], makers[n][2]().fit(L.iloc[a], yL[a]).predict(L.iloc[b])) for a, b in fl]))
        leak.append(row)
        print(f"  leakage check {n:28s} grouped {row['grouped']:.4f}  random {row['random']:.4f}", flush=True)

    # ---- 4. repeated grouped splits of training + validation (same settings, nothing retuned)
    rep_names = list(dict.fromkeys(["Mean", "Mutant residue mean", REF, "Ridge (F1)", "Ridge (F2)", selected]))
    repeats = {n: [] for n in rep_names}
    gpool = pool[c["group"]].to_numpy()
    for r in range(cfg["n_repeats"]):
        a, b = next(GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=1000 + r).split(pool, groups=gpool))
        for n in rep_names:
            repeats[n].append(mae(ypool[b], makers[n][2]().fit(pool.iloc[a], ypool[a]).predict(pool.iloc[b])))
        print(f"  repeat {r + 1}/{cfg['n_repeats']}  {time.time() - t0:.0f}s", flush=True)
    rep_out = {}
    for n in rep_names:
        v = np.array(repeats[n])
        rep_out[n] = {"mae": [round(x, 4) for x in v], "mae_mean": float(v.mean()),
                      "n_better_mean": int((v < np.array(repeats["Mean"])).sum()),
                      "n_better_ref": int((v < np.array(repeats[REF])).sum())}

    # ---- 5. error analysis of the selected model on the test set
    p, e = pte[selected], pte[selected] - yte
    origin = np.where(te[c["group"]].str.fullmatch(r"\d+"), "natural", "designed")
    stype = np.array(["to " + CLASS.get(a, "other") for a in te[c["mutated"]]])
    size = pd.cut(np.abs(yte), [-0.01, 0.5, 1, 2, np.inf], labels=["|ddG| under 0.5", "0.5 to 1", "1 to 2", "over 2"]).astype(str)

    def slices(labels, order=None):
        out = []
        for lab in (order or sorted(set(labels))):
            if not (labels == lab).any():
                continue
            i = labels == lab
            out.append({"slice": str(lab), "n": int(i.sum()), "mae": mae(yte[i], p[i]), "bias": float(e[i].mean()),
                        "mae_mean": mae(yte[i], pte["Mean"][i]), "mae_ref": mae(yte[i], pte[REF][i])})
        return out
    by_cluster = sorted(slices(te[c["group"]].to_numpy()), key=lambda s: -s["n"])
    for s in by_cluster:
        s["origin"] = "natural" if s["slice"].isdigit() else "designed"
        s["proteins"] = int(te.loc[te[c["group"]] == s["slice"], c["protein"]].nunique())
    residuals = {"by_type": slices(stype), "by_origin": slices(origin), "by_size": slices(np.asarray(size), ["|ddG| under 0.5", "0.5 to 1", "1 to 2", "over 2"]), "by_cluster": by_cluster,
                 "pred_sd": float(p.std()), "y_sd": float(yte.std())}

    # ---- 6. permutation importance of the selected model, on validation rows
    fit_sel = makers[selected][2]().fit(tr, ytr)
    vi = np.sort(rng.choice(len(va), min(cfg["importance_rows"], len(va)), replace=False))
    V, yV = va.iloc[vi].reset_index(drop=True), yva[vi]
    base = mae(yV, fit_sel.predict(V))
    importance = []
    for col in FS[sel_fs]["cat"] + FS[sel_fs]["num"]:
        inc = []
        for r in range(cfg["importance_repeats"]):
            Vp = V.copy()
            Vp[col] = Vp[col].sample(frac=1, random_state=seed + r).to_numpy()
            inc.append(mae(yV, fit_sel.predict(Vp)) - base)
        importance.append({"feature": col, "mae_increase": float(np.mean(inc)), "sd": float(np.std(inc))})
    importance.sort(key=lambda d: -d["mae_increase"])
    print(f"  importance done  {time.time() - t0:.0f}s", flush=True)

    # ---- 7. rows for the website: test predictions, sandbox rows, parity checks
    shown = np.sort(rng.choice(len(te), min(cfg["test_rows_shown"], len(te)), replace=False))
    show_models = ["Mean", REF] + learned
    test_rows = [{"mut": te[c["mutation"]].iat[j], "wt": te[c["protein"]].iat[j], "group": te[c["group"]].iat[j],
                  "stratum": stype[j], "designed": int(origin[j] == "designed"), "y": round(float(yte[j]), 3),
                  "pred": {n: round(float(pte[n][j]), 3) for n in show_models}} for j in shown]
    pd.DataFrame({"protein": te[c["protein"]], "cluster": te[c["group"]], "mutation": te[c["mutation"]], "ddG": yte,
                  **{n: np.round(pte[n], 4) for n in makers}}).to_csv(ROOT / "results" / "test_predictions.csv.gz", index=False, compression="gzip")
    sandbox = {"rows": [], "x_check": [], "check": None}
    if c.get("sequence") and c["sequence"] in tr.columns:
        allp = pd.concat([tr.assign(_test=False), va.assign(_test=False), te.assign(_test=True)], ignore_index=True)
        si = np.sort(np.random.default_rng(seed + 1).choice(len(allp), min(cfg["sandbox_rows"], len(allp)), replace=False))
        S = allp.iloc[si].reset_index(drop=True)
        sandbox["rows"] = [{"seq": S[c["sequence"]].iat[i], "mut": S[c["mutation"]].iat[i], "y": round(float(S[T].iat[i]), 4),
                            "group": S[c["group"]].iat[i], "test": bool(S["_test"].iat[i])} for i in range(len(S))]
        sandbox["x_check"] = [[float(v) for v in row] for row in design_f1(S.iloc[:40], c)]
        a, b = ~S["_test"].to_numpy(), S["_test"].to_numpy()
        ys = S[T].to_numpy(float)
        sm = pipeline(FS["F1"], Ridge(alpha=cfg["sandbox_alpha"])).fit(S[a], ys[a])
        sandbox["check"] = {"alpha": cfg["sandbox_alpha"], "test_mae": mae(ys[b], sm.predict(S[b]))}

    allrows = pd.concat([tr, va, te], ignore_index=True)
    yall = allrows[T].to_numpy(float)
    des = ~allrows[c["group"]].str.fullmatch(r"\d+").to_numpy()
    lo, step = float(np.floor(yall.min())), 0.25
    edges = np.arange(lo, np.ceil(yall.max()) + step, step)
    cl = allrows[c["group"]].value_counts()
    facts = cfg.get("facts", {})
    summary = {
        **facts, "n_kept": int(len(allrows)), "n_proteins": int(allrows[c["protein"]].nunique()), "n_clusters": int(len(cl)),
        "largest_cluster": str(cl.index[0]), "largest_cluster_n": int(cl.iloc[0]), "n_designed": int(des.sum()),
        "n_train": len(tr), "n_val": len(va), "n_test": len(te),
        "c_train": int(tr[c["group"]].nunique()), "c_val": int(va[c["group"]].nunique()), "c_test": int(te[c["group"]].nunique()),
        "target_mean": float(yall.mean()), "target_sd": float(yall.std()), "target_min": float(yall.min()), "target_max": float(yall.max()),
        "share_small": float((np.abs(yall) < 0.5).mean()), "share_destabilizing": float((yall < 0).mean()),
        "len_min": int(allrows[c["length"]].min()), "len_max": int(allrows[c["length"]].max()),
        "hist": {"lo": lo, "step": step, "designed": np.histogram(yall[des], edges)[0].tolist(), "natural": np.histogram(yall[~des], edges)[0].tolist()},
        "features": {k: v["cat"] + v["num"] for k, v in FS.items()},
    }
    import sklearn, scipy
    out = {"config": cfg, "generated": time.strftime("%Y-%m-%d"), "selected": selected, "ref": REF,
           "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__,
                        "scikit-learn": sklearn.__version__, "scipy": scipy.__version__},
           "summary": summary, "models": models, "leak": leak, "repeats": rep_out, "residuals": residuals,
           "importance": importance, "test_rows": test_rows, "sandbox": sandbox}
    (ROOT / "results" / "evaluation.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    (ROOT / "docs" / "results.js").write_text("window.RESULTS = " + json.dumps(out, separators=(",", ":"), ensure_ascii=False) + ";\n", encoding="utf-8")

    print(f"\n{'model':32s} {'val MAE':>8s} {'test MAE':>9s} {'95% interval':>17s} {'vs mean':>8s} {'vs table [95% CI]':>24s}")
    for m in models:
        print(f"{m['name']:32s} {m['val_mae']:8.4f} {m['test_mae']:9.4f} [{m['test_mae_ci'][0]:.3f}, {m['test_mae_ci'][1]:.3f}] "
              f"{m['diff_vs_mean']:+8.4f} {m['diff_vs_ref']:+8.4f} [{m['diff_vs_ref_ci'][0]:+.3f}, {m['diff_vs_ref_ci'][1]:+.3f}]")
    print(f"selected: {selected} | done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "docs").mkdir(exist_ok=True)
    main()
