import re, math
from pathlib import Path
import pandas as pd

def discover_csvs(folder: Path) -> dict:
    files = {}
    if folder.is_dir():
        for p in folder.glob("*.csv"):
            files[p.stem] = str(p.resolve())
    return files

def _cond_col(df):
    for c in df.columns:
        if c.lower() in ("group","condition","treatment","label","cond"): return c
    for c in df.columns:
        if re.search("cond", c, re.I): return c
    return None

def _num_col(df):
    for c in df.columns:
        if df[c].dtype.kind in "fi": return c
    return None

def parse_pe_3b(path):
    """
    Parse Pe Fig 3B data (timecourse per condition) and estimate delta_r/day.
    Uses Mean_FoldChange column when present.
    """
    import pandas as pd, numpy as np, re, math
    from pathlib import Path

    p = Path(path)
    if not p.exists():
        print(f"[WARN] 3BPe not found: {p}")
        return {}

    df = pd.read_csv(p)

    # --- column detection ---
    cond = None
    time_col = None
    val_col = None
    for c in df.columns:
        if re.search(r"(cond|group|label|treat|type)", c, re.I):
            cond = c
        elif re.search(r"(time|hour|h)", c, re.I):
            time_col = c
        elif re.search(r"(mean|value|fold)", c, re.I):      # prefer Mean_FoldChange
            val_col = c
    if not all([cond, time_col, val_col]):
        print(f"[WARN] 3BPe columns not recognized: {list(df.columns)}")
        return {}

    # --- convert hours to days ---
    df["_t_days"] = df[time_col] / 24.0 if df[time_col].max() > 10 else df[time_col]

    # --- select Control and TAM only ---
    conds = df[cond].astype(str)
    ctrl = df[conds.str.contains(r"ctrl|control|vehicle|isotype", case=False, regex=True)]
    tam  = df[conds.str.contains(r"TAM", case=False, regex=True)]
    if ctrl.empty or tam.empty:
        print(f"[WARN] 3BPe missing control or TAM rows.")
        return {}

    # --- fit exponential growth r = slope(log(C)) vs t ---
    def fit_rate(sub):
        t = sub["_t_days"].values
        y = np.log(np.clip(sub[val_col].values, 1e-9, None))
        if len(t) < 2:
            return np.nan
        A = np.vstack([t, np.ones_like(t)]).T
        slope, _ = np.linalg.lstsq(A, y, rcond=None)[0]
        return slope

    r_ctrl = fit_rate(ctrl)
    r_tam  = fit_rate(tam)
    if np.isnan(r_ctrl) or np.isnan(r_tam):
        print("[WARN] could not fit exponential to 3BPe data")
        return {}

    delta_r = float(r_tam - r_ctrl)
    fold = math.exp(delta_r * 3.0)  # expected 72 h fold-change difference

    print(f"[INFO] r_ctrl = {r_ctrl:.3f}  r_tam = {r_tam:.3f}  Δr = {delta_r:.3f}")
    return {"fold_TAM": fold, "delta_r_TAM_per_day": delta_r}


def parse_pe_2d(path):
    p = Path(path)
    if not p.exists(): return {"TAM": {}, "M2": {}}
    df = pd.read_csv(p)
    cond = _cond_col(df)
    out = {"TAM": {}, "M2": {}}
    for col in df.columns:
        if re.search("IL|TNF", col, re.I):
            out["TAM"][col] = df[df[cond].str.contains("TAM", case=False, na=False)][col].mean()
            out["M2"][col]  = df[df[cond].str.contains("M2", case=False, na=False)][col].mean()
    return out

def parse_qian_sigma_ccl2(path):
    """
    Extract anti-CCL2/control ratio from Qian Fig1X CSVs.
    Matches flexible label patterns like 'Isotype', 'Ctrl', 'anti-CCL2 Ab', etc.
    """
    from pathlib import Path
    import pandas as pd, re

    p = Path(path)
    if not p.exists():
        print(f"[WARN] Qian CSV not found: {p}")
        return None

    df = pd.read_csv(p)
    cond = None
    for c in df.columns:
        if re.search(r"(cond|group|label|type|treat)", c, re.I):
            cond = c; break
    val = None
    for c in df.columns:
        if df[c].dtype.kind in "fi":
            val = c; break
    if cond is None or val is None:
        print(f"[WARN] Qian columns not recognized: {list(df.columns)}")
        return None

    conds = df[cond].astype(str)
    ctrl = df[conds.str.contains(r"ctrl|control|vehicle|isotype", case=False, regex=True)][val]
    anti = df[conds.str.contains(r"anti[-_ ]?CCL2|CCL2[-_ ]?Ab|αCCL2", case=False, regex=True)][val]

    if ctrl.empty or anti.empty:
        print(f"[WARN] Could not find control/anti rows in {p.name}")
        return None

    ctrl_mean = float(ctrl.mean())
    anti_mean = float(anti.mean())
    if ctrl_mean <= 0: return None

    ratio = anti_mean / ctrl_mean
    return ratio

def parse_qian_kappa_vegf(path):
    p = Path(path)
    if not p.exists(): return None
    df = pd.read_csv(p)
    cond, val = _cond_col(df), _num_col(df)
    wt = df[df[cond].str.contains("WT", case=False, na=False)][val].mean()
    ko = df[df[cond].str.contains("KO", case=False, na=False)][val].mean()
    return float(ko/wt) if wt else None
