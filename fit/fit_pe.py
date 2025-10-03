
# fit/fit_pe.py — minimal fitter for Pe proliferation data
import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model.equations import simulate_proliferation
from fit.utils import wls_loss, sem_to_weight

# Data loader
def load_pe(csv_path):
    df = pd.read_csv(csv_path)
    # Accept either "Time" or "Time_h"
    if "Time_h" not in df.columns and "Time" in df.columns:
        df = df.rename(columns={"Time": "Time_h"})
    # If SEM is split as SEM+ / SEM-, collapse to one SEM
    if "SEM" not in df.columns:
        if "SEM+" in df.columns and "SEM-" in df.columns:
            df["SEM"] = 0.5*(df["SEM+"].abs() + df["SEM-"].abs())
        else:
            df["SEM"] = np.nan
    return df


def compute_relfold(df):
    # Pivot by time and condition
    # RelFold(cond,t) = Mean(cond,t)/Mean(Control,t)
    times = sorted(df["Time_h"].unique())
    conds = sorted(df["Condition"].unique())
    # build dict time -> {cond: (mean, sem)}
    out_rows = []
    for t in times:
        sub = df[df["Time_h"]==t]
        ctrl = sub[sub["Condition"].str.lower()=="control"]
        if ctrl.empty:
            continue
        ctrl_mean = float(ctrl["Mean"].values[0])
        ctrl_sem = float(ctrl["SEM"].values[0]) if not np.isnan(ctrl["SEM"].values[0]) else np.nan
        for cond in conds:
            row = sub[sub["Condition"]==cond]
            if row.empty: 
                continue
            m = float(row["Mean"].values[0])
            s = float(row["SEM"].values[0]) if not np.isnan(row["SEM"].values[0]) else np.nan
            rel = m/ctrl_mean if ctrl_mean!=0 else np.nan
            # propagate sem approximately: rel * sqrt( (s/m)^2 + (ctrl_sem/ctrl_mean)^2 )
            if not np.isnan(s) and m>0 and ctrl_mean>0 and not np.isnan(ctrl_sem):
                rel_sem = rel*np.sqrt((s/max(m,1e-9))**2 + (ctrl_sem/max(ctrl_mean,1e-9))**2)
            else:
                rel_sem = np.nan
            out_rows.append({"Time_h":t,"Condition":cond,"RelFold":rel,"RelSEM":rel_sem})
    return pd.DataFrame(out_rows)

def fit_params(pe_path, results_dir):
    df = load_pe(pe_path)
    rel = compute_relfold(df)
    # Times to simulate (ensures sorted unique)
    times = np.array(sorted(rel["Time_h"].unique()), dtype=float)
    # conditions
    conds = [c for c in sorted(rel["Condition"].unique()) if c.lower()!="control"]
    # Observations and weights
    obs = {cond: rel[rel["Condition"]==cond].sort_values("Time_h")["RelFold"].values for cond in conds}
    sem = {cond: rel[rel["Condition"]==cond].sort_values("Time_h")["RelSEM"].values for cond in conds}
    wts = {cond: np.where(np.isnan(sem[cond]), 1.0, 1.0/(sem[cond]**2 + 1e-9)) for cond in conds}

    for cond in conds:
        w = wts[cond]
        wts[cond] = np.maximum(w, 1.0)  # floor so huge SEM doesn't zero-out the fit

    # Simple parameterization:
    # params: r0, gamma2, K2, s_M2, s_TAM, theta, (optional gamma1,K1= small)
    params = dict(r0=np.log(2)/28.0, gamma2=0.5, K2=1.0, gamma1=0.2, K1=1.0)
    scales = dict(M2=1.0, TAM=1.0, theta=0.3)

    # Define loss function
    def loss(p):
        r0, gamma2, K2, s_M2, s_TAM, theta = p
        pars = dict(r0=r0, gamma2=gamma2, K2=K2, gamma1=params["gamma1"], K1=params["K1"])
        # simulate control
        C_ctrl = simulate_proliferation(times, pars, dict(S_M1=0.0,S_M2=0.0,C0=1.0))
        total = 0.0
        for cond in conds:
            if cond.lower()=="m2":
                scen = dict(S_M1=0.0, S_M2=s_M2, C0=1.0)
            else: # TAM
                scen = dict(S_M1=theta*s_TAM, S_M2=(1.0-theta)*s_TAM, C0=1.0)
            C_cond = simulate_proliferation(times, pars, scen)
            pred = C_cond / C_ctrl
            total += np.sum(wts[cond]*(pred - obs[cond])**2)
        return total

    # Crude grid search + local tweaks (avoid SciPy dependency)
    best = None
    best_val = 1e99
    # ----- search ranges (wider) -----
    r0_range = np.linspace(np.log(2) / 40, np.log(2) / 18, 12)
    gamma2_range = np.linspace(0.0, 3.0, 13)
    K2_range = np.linspace(0.05, 10.0, 12)
    sM2_range = np.linspace(0.1, 4.0, 12)
    sTAM_range = np.linspace(0.1, 4.0, 12)
    theta_range = np.linspace(0.1, 0.7, 7)

    # ----- weight floor so huge SEM can't zero-out points -----
    for cond in list(wts.keys()):
        wts[cond] = np.maximum(wts[cond], 1.0)

    # ----- grid search -----
    best = None
    best_val = 1e99
    for r0 in r0_range:
        for g2 in gamma2_range:
            for k2 in K2_range:
                for sM2 in sM2_range:
                    for sT in sTAM_range:
                        for th in theta_range:
                            val = loss([r0, g2, k2, sM2, sT, th])
                            if val < best_val:
                                best_val = val
                                best = [r0, g2, k2, sM2, sT, th]

    # Prepare outputs
    r0,g2,k2,sM2,sT,th = best
    best_params = dict(r0=float(r0), gamma2=float(g2), K2=float(k2), gamma1=float(params["gamma1"]), K1=float(params["K1"]),
                       s_M2=float(sM2), s_TAM=float(sT), theta=float(th))
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir,"pe_fit.json"),"w") as f:
        json.dump(best_params, f, indent=2)

    # Build plot
    pars = dict(r0=best_params["r0"], gamma2=best_params["gamma2"], K2=best_params["K2"],
                gamma1=best_params["gamma1"], K1=best_params["K1"])
    C_ctrl = simulate_proliferation(times, pars, dict(S_M1=0.0,S_M2=0.0,C0=1.0))
    fig, ax = plt.subplots(figsize=(6,4))
    # plot data points
    for cond, color, marker in zip(conds, ["tab:green","tab:orange"], ["o","s"]):
        y = obs[cond]
        ysem = sem[cond]
        ax.errorbar(times, y, yerr=ysem, fmt=marker, label=f"Data {cond}", capsize=3)
        if cond.lower()=="m2":
            scen = dict(S_M1=0.0, S_M2=best_params["s_M2"], C0=1.0)
        else:
            scen = dict(S_M1=best_params["theta"]*best_params["s_TAM"],
                        S_M2=(1.0-best_params["theta"])*best_params["s_TAM"], C0=1.0)
        C_cond = simulate_proliferation(times, pars, scen)
        pred = C_cond / C_ctrl
        ax.plot(times, pred, '-', label=f"Model {cond}")
    ax.set_title("Pe 2022 — Rel. fold-change vs Control (WST)")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Rel. Fold vs Control")
    ax.legend()
    plt.tight_layout()
    figpath = os.path.join(results_dir, "pe_fit.png")
    plt.savefig(figpath, dpi=160)
    plt.close(fig)

    return best_params, figpath

if __name__ == "__main__":
    # Resolve data path
    try:
        from config import get_pe_path
        pe_path = get_pe_path()
    except Exception:
        pe_path = os.path.join(os.path.dirname(__file__), "..", "data", "Pe_Proliferation.csv")
    results_dir = os.path.join(os.path.dirname(__file__), "..", "results", "figures")
    params, figpath = fit_params(pe_path, results_dir)
    print("Saved:", figpath)
    print("Params:", json.dumps(params, indent=2))
