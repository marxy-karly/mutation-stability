"""Extra features for evaluate.py. Nothing here is fitted to the data, so nothing can leak.

F3 = F2 + sequence-derived numbers, available for every protein:
    how the swap changes hydrophobicity, size, charge and helix propensity;
    where the residue sits (relative position, distance to the nearest end);
    how hydrophobic its neighbourhood is, overall and on the residues that would face the
    same way in a helix (i+-3, i+-4) or a strand (i+-2). A hydrophobic face suggests a buried residue.

F4 = F3 + burial in the AlphaFold model of the unmutated protein:
    number of residues within 8 and 10 Angstrom, half-sphere exposure, AlphaFold confidence (pLDDT),
    a helix/strand geometry measure, and the residue's burial rank inside its protein.
    Proteins without a structure get missing values, filled with the training median inside the model
    pipeline, plus a 0/1 column st_found.

Structures: AlphaFold_model_PDBs.zip from Zenodo record 7401275, unzipped into data/structures/
(any folder layout). Structures are matched to proteins by SEQUENCE, not file name: the
wild-type sequence is rebuilt from aa_seq by putting the original residue back at the mutated position.

Check coverage without running the models:
    python src/features.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
AAS = "ACDEFGHIKLMNPQRSTVWY"
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
         "MSE": "M", "HSD": "H", "HSE": "H", "HIE": "H", "HID": "H"}
# Kyte & Doolittle (1982) hydrophobicity
HYD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6,
           S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
# Zamyatnin (1972) residue volume, cubic Angstrom
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1, H=153.2, I=166.7, L=166.7, K=168.6, M=162.9,
           F=189.9, P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)
CHG = dict(K=1.0, R=1.0, D=-1.0, E=-1.0)
# Pace & Scholtz (1998) helix propensity, kcal/mol relative to alanine (higher = worse in a helix)
HELIX = dict(A=0.0, L=0.21, R=0.21, M=0.24, K=0.26, Q=0.39, E=0.40, I=0.41, W=0.49, S=0.50, Y=0.53, F=0.54, H=0.61, V=0.61,
             N=0.65, T=0.66, C=0.68, D=0.69, G=1.00, P=3.16)

SEQ_COLS = ["hyd_wt", "hyd_mut", "hyd_delta", "vol_wt", "vol_mut", "vol_delta", "charge_delta", "helix_delta",
            "rel_position", "dist_to_end", "win_hyd_4", "win_hyd_8", "win_charge_4", "helix_face_hyd", "strand_face_hyd"]
# what each column tells the model, used to group the permutation importance
SWAP_SEQ_COLS = ["hyd_wt", "hyd_mut", "hyd_delta", "vol_wt", "vol_mut", "vol_delta", "charge_delta", "helix_delta"]
POS_SEQ_COLS = ["rel_position", "dist_to_end"]
NEIGH_SEQ_COLS = ["win_hyd_4", "win_hyd_8", "win_charge_4", "helix_face_hyd", "strand_face_hyd"]
STRUCT_COLS = ["st_found", "st_contacts_8", "st_contacts_10", "st_hse_up", "st_plddt", "st_local_span", "st_burial_rank",
               "st_rsa", "st_neigh_hyd", "st_centroid_dist", "x_hyd_burial", "x_vol_burial", "x_helix_burial"]
# Tien et al. (2013) theoretical maximum accessible surface area per residue, square Angstrom
MAX_ASA = dict(A=129.0, R=274.0, N=195.0, D=193.0, C=167.0, Q=225.0, E=223.0, G=104.0, H=224.0, I=197.0, L=201.0, K=236.0,
               M=224.0, F=240.0, P=159.0, S=155.0, T=172.0, W=285.0, Y=263.0, V=174.0)
RADIUS = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}   # heavy-atom radii for the surface calculation
PROBE = 1.4
STRUCT_DESC = {
    "st_found": "1 if the protein has a matched AlphaFold model",
    "st_contacts_8": "residues with a side-chain anchor (C-beta, C-alpha for glycine) within 8 Å",
    "st_contacts_10": "the same within 10 Å",
    "st_hse_up": "C-alpha atoms within 13 Å on the side the side chain points to (half-sphere exposure)",
    "st_plddt": "AlphaFold's confidence at that residue",
    "st_local_span": "distance from residue i−2 to i+2: about 6 Å in a helix, 12 Å or more in a strand",
    "st_burial_rank": "rank of st_contacts_10 within the protein, 0 = most exposed, 1 = most buried",
    "st_rsa": "relative solvent-accessible surface area (0 = fully buried, 1 = fully exposed), Shrake-Rupley",
    "st_neigh_hyd": "mean hydrophobicity of the residues within 10 Å in 3D (not in sequence)",
    "st_centroid_dist": "distance from the protein's centre, divided by its radius of gyration",
    "x_hyd_burial": "hydrophobicity change × burial (1 − relative surface area)",
    "x_vol_burial": "volume change × burial",
    "x_helix_burial": "helix-propensity change × burial",
}
INTERACTIONS = {"x_hyd_burial": "hyd_delta", "x_vol_burial": "vol_delta", "x_helix_burial": "helix_delta"}


def _scale(seq, table, default=0.0):
    return np.array([table.get(a, default) for a in seq], dtype=float)


def residue_table(wt):
    """Per-residue sequence features of one wild-type sequence: window averages that do not include the residue itself."""
    L, h, q = len(wt), _scale(wt, HYD), _scale(wt, CHG)
    out = {}
    for name, vec, w, agg in (("win_hyd_4", h, 4, np.mean), ("win_hyd_8", h, 8, np.mean), ("win_charge_4", q, 4, np.sum)):
        vals = np.empty(L)
        for i in range(L):
            idx = [j for j in range(max(0, i - w), min(L, i + w + 1)) if j != i]
            vals[i] = agg(vec[idx]) if idx else 0.0
        out[name] = vals   # mean hydrophobicity, or net charge, of the neighbours
    for name, ks in (("helix_face_hyd", (3, 4)), ("strand_face_hyd", (2,))):
        vals = np.empty(L)
        for i in range(L):
            idx = [i + s * k for k in ks for s in (-1, 1) if 0 <= i + s * k < L]
            vals[i] = h[idx].mean() if idx else 0.0
        out[name] = vals
    return pd.DataFrame(out)


# ----------------------------------------------------------------------------- structures
def read_pdb(path):
    """First model, first chain. Returns sequence, C-alpha and side-chain anchor coordinates, and pLDDT (B-factor) per residue."""
    res, order, chain = {}, [], None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM"):
                continue
            ch = line[21]
            if chain is None:
                chain = ch
            if ch != chain:
                continue
            key = (int(line[22:26]), line[26])
            if key not in res:
                aa = THREE.get(line[17:20].strip())
                if aa is None:
                    continue
                res[key] = {"aa": aa, "b": float(line[60:66] or 0), "heavy": []}
                order.append(key)
            atom = line[12:16].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if atom in ("CA", "CB"):
                res[key][atom] = xyz
            el = (line[76:78].strip() or atom[0]).upper()
            if el in RADIUS:
                res[key]["heavy"].append((xyz, RADIUS[el]))
    order = [k for k in order if "CA" in res[k]]
    if not order:
        return None
    seq = "".join(res[k]["aa"] for k in order)
    ca = np.array([res[k]["CA"] for k in order])
    cb = np.array([res[k].get("CB", res[k]["CA"]) for k in order])
    heavy = [(i, a[0], a[1]) for i, k in enumerate(order) for a in res[k]["heavy"]]
    return {"seq": seq, "ca": ca, "cb": cb, "plddt": np.array([res[k]["b"] for k in order]),
            "atom_res": np.array([h[0] for h in heavy], dtype=int), "atom_xyz": np.array([h[1] for h in heavy], dtype=float),
            "atom_r": np.array([h[2] for h in heavy], dtype=float)}


def _sphere(n=96):
    """Evenly spread points on a unit sphere (golden-spiral)."""
    i = np.arange(n) + 0.5
    phi, theta = np.arccos(1 - 2 * i / n), np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])


SPHERE = _sphere()


def residue_sasa(s):
    """Shrake-Rupley solvent-accessible surface area per residue, heavy atoms only."""
    xyz, r, res = s["atom_xyz"], s["atom_r"] + PROBE, s["atom_res"]
    L = len(s["seq"])
    if len(xyz) == 0:
        return np.full(L, np.nan)
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    area = np.zeros(len(xyz))
    for i in range(len(xyz)):
        nb = np.where((d[i] < r[i] + r) & (np.arange(len(xyz)) != i))[0]
        pts = xyz[i] + r[i] * SPHERE
        if len(nb):
            free = (np.linalg.norm(pts[:, None, :] - xyz[nb][None, :, :], axis=-1) >= r[nb][None, :]).all(1)
        else:
            free = np.ones(len(SPHERE), bool)
        area[i] = 4 * np.pi * r[i] ** 2 * free.mean()
    return np.bincount(res, weights=area, minlength=L)


def structure_table(s):
    """Per-residue burial features of one structure."""
    ca, cb, L = s["ca"], s["cb"], len(s["seq"])
    d = np.linalg.norm(cb[:, None, :] - cb[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    c8, c10 = (d < 8).sum(1).astype(float), (d < 10).sum(1).astype(float)
    # half-sphere exposure: C-alpha atoms within 13 A in the half-space the C-alpha -> C-beta vector points to
    v = cb - ca
    for i in np.where(np.linalg.norm(v, axis=1) < 1e-6)[0]:   # glycine: point away from the two neighbouring C-alphas
        if 0 < i < L - 1:
            v[i] = ca[i] - (ca[i - 1] + ca[i + 1]) / 2
    norm = np.linalg.norm(v, axis=1)
    hse = np.full(L, np.nan)
    dca = ca[None, :, :] - ca[:, None, :]
    within = np.linalg.norm(dca, axis=-1) < 13
    np.fill_diagonal(within, False)
    ok = norm > 1e-6
    hse[ok] = ((np.einsum("ijk,ik->ij", dca[ok], v[ok]) > 0) & within[ok]).sum(1)
    span = np.full(L, np.nan)
    if L >= 5:
        span[2:L - 2] = np.linalg.norm(ca[4:] - ca[:-4], axis=1)
    rank = pd.Series(c10).rank(pct=True).to_numpy()
    rank = (rank - rank.min()) / (rank.max() - rank.min()) if rank.max() > rank.min() else np.full(L, 0.5)
    rsa = np.clip(residue_sasa(s) / np.array([MAX_ASA[a] for a in s["seq"]]), 0, 1)
    h = _scale(s["seq"], HYD)
    near = d < 10
    neigh_hyd = np.array([h[near[i]].mean() if near[i].any() else 0.0 for i in range(L)])
    cen = ca.mean(0)
    rg = np.sqrt(((ca - cen) ** 2).sum(1).mean()) or 1.0
    cdist = np.linalg.norm(cb - cen, axis=1) / rg
    # the interaction columns are filled per row in add_features (they need the swap); NaN placeholders here
    nan = np.full(L, np.nan)
    return pd.DataFrame({"st_found": 1.0, "st_contacts_8": c8, "st_contacts_10": c10, "st_hse_up": hse, "st_plddt": s["plddt"],
                         "st_local_span": span, "st_burial_rank": rank, "st_rsa": rsa, "st_neigh_hyd": neigh_hyd,
                         "st_centroid_dist": cdist, "x_hyd_burial": nan, "x_vol_burial": nan, "x_helix_burial": nan})


def load_structures(folder):
    files = sorted(p for p in Path(folder).rglob("*") if p.is_file() and p.suffix.lower() in (".pdb", ".ent"))
    out = []
    for p in files:
        s = read_pdb(p)
        if s:
            s["file"] = p.name
            out.append(s)
    return out


def match(wt, structs, by_seq, by_len, max_mismatch):
    """Find the structure for one wild-type sequence. Returns (structure, offset, how) or (None, 0, 'none').
    offset maps a 0-based sequence index i to structure index i + offset."""
    if wt in by_seq:
        return by_seq[wt], 0, "exact"
    best, bm = None, max_mismatch + 1
    a = np.frombuffer(wt.encode(), dtype=np.uint8)
    for s in by_len.get(len(wt), []):
        mm = int((np.frombuffer(s["seq"].encode(), dtype=np.uint8) != a).sum())
        if mm < bm:
            best, bm = s, mm
    if best is not None:
        return best, 0, "near"
    for s in structs:   # one contains the other: extra tags or linkers on one side
        if len(s["seq"]) >= 20 and len(wt) >= 20:
            if wt in s["seq"]:
                return s, s["seq"].find(wt), "offset"
            if s["seq"] in wt:
                return s, -wt.find(s["seq"]), "offset"
    return None, 0, "none"


# ----------------------------------------------------------------------------- main entry
def add_features(parts, cfg, log=print):
    """Adds SEQ_COLS (always) and STRUCT_COLS (when structures are found) to every split, in place.
    Returns a dict describing what was added, for the website."""
    c = cfg["columns"]
    S, P, O, Mu = c["sequence"], c["position"], c["original"], c["mutated"]
    info = {"seq_cols": SEQ_COLS, "struct_cols": [], "struct_desc": STRUCT_DESC}
    folder = ROOT / cfg.get("structures", "data/structures")
    structs = load_structures(folder) if folder.exists() else []
    by_seq, by_len = {}, {}
    for s in structs:
        by_seq.setdefault(s["seq"], s)
        by_len.setdefault(len(s["seq"]), []).append(s)
    max_mm = int(cfg.get("structure_max_mismatch", 3))
    cache, how_count, struct_cache = {}, {}, {}
    seq_bad = 0
    for part, df in parts.items():
        pos = pd.to_numeric(df[P], errors="coerce").astype(int).to_numpy()
        seqs, orig, new = df[S].astype(str).to_numpy(), df[O].to_numpy(), df[Mu].to_numpy()
        wt = np.empty(len(df), dtype=object)
        for i in range(len(df)):
            q, p = seqs[i], pos[i] - 1
            if 0 <= p < len(q) and q[p] == new[i]:
                wt[i] = q[:p] + orig[i] + q[p + 1:]
            else:
                wt[i] = q
                seq_bad += 1
        L = np.array([len(x) for x in wt], dtype=float)
        df["hyd_wt"], df["hyd_mut"] = df[O].map(HYD).astype(float), df[Mu].map(HYD).astype(float)
        df["hyd_delta"] = df["hyd_mut"] - df["hyd_wt"]
        df["vol_wt"], df["vol_mut"] = df[O].map(VOL).astype(float), df[Mu].map(VOL).astype(float)
        df["vol_delta"] = df["vol_mut"] - df["vol_wt"]
        df["charge_delta"] = df[Mu].map(CHG).fillna(0.0) - df[O].map(CHG).fillna(0.0)
        df["helix_delta"] = df[Mu].map(HELIX).astype(float) - df[O].map(HELIX).astype(float)
        df["rel_position"] = (pos - 1) / np.maximum(L - 1, 1)
        df["dist_to_end"] = np.minimum(pos - 1, L - pos).astype(float)
        win = np.full((len(df), 5), np.nan)
        st = np.full((len(df), len(STRUCT_COLS)), np.nan)
        st[:, 0] = 0.0
        groups = pd.Series(np.arange(len(df))).groupby(wt).indices
        for w, idx in groups.items():
            if w not in cache:
                cache[w] = residue_table(w).to_numpy()
            rows = pos[idx] - 1
            okr = (rows >= 0) & (rows < len(w))
            win[idx[okr]] = cache[w][rows[okr]]
            if structs:
                if w not in struct_cache:
                    s, off, how = match(w, structs, by_seq, by_len, max_mm)
                    struct_cache[w] = (structure_table(s).to_numpy() if s else None, off, how, s["file"] if s else None)
                    how_count[how] = how_count.get(how, 0) + 1
                tab, off, how, _ = struct_cache[w]
                if tab is not None:
                    j = rows + off
                    okj = (j >= 0) & (j < len(tab))
                    st[idx[okj]] = tab[j[okj]]
        for k, col in enumerate(["win_hyd_4", "win_hyd_8", "win_charge_4", "helix_face_hyd", "strand_face_hyd"]):
            df[col] = win[:, k]
        if structs:
            for k, col in enumerate(STRUCT_COLS):
                df[col] = st[:, k]
            for col, src in INTERACTIONS.items():
                df[col] = df[src] * (1 - df["st_rsa"])
    if seq_bad:
        log(f"  note: {seq_bad} rows did not have the new residue at the stated position; their own sequence was used as the wild type")
    if structs:
        info["struct_cols"] = STRUCT_COLS
        cov = {}
        for part, df in parts.items():
            found = df["st_found"].to_numpy() == 1
            cov[part] = {"rows": int(found.sum()), "rows_share": float(found.mean()),
                         "proteins": int(df.loc[found, c["protein"]].nunique()), "proteins_total": int(df[c["protein"]].nunique())}
        info["coverage"] = cov
        info["structures_read"] = len(structs)
        info["match_counts"] = how_count
        log(f"  structures: read {len(structs)} files from {folder.relative_to(ROOT)}; wild types matched {how_count}")
        for part, v in cov.items():
            log(f"    {part:10s} {v['rows_share']:6.1%} of rows, {v['proteins']}/{v['proteins_total']} proteins have a structure")
        if sum(v["rows"] for v in cov.values()) == 0:
            info["struct_cols"] = []
            log("  no row matched a structure, so F4 is skipped")
    else:
        log(f"  no structures in {folder.relative_to(ROOT) if folder.is_relative_to(ROOT) else folder}: F4 is skipped "
            "(download AlphaFold_model_PDBs.zip from Zenodo 7401275 and unzip it there)")
    return info


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from evaluate import load
    cfg = json.loads((ROOT / "config.json").read_text(encoding="utf-8"))
    parts = load(cfg)
    info = add_features(parts, cfg)
    print(json.dumps({k: v for k, v in info.items() if k not in ("struct_desc",)}, indent=1))
