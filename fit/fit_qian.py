
# fit/fit_qian.py — minimal fitter for Qian recruitment ratios
import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model.equations import simulate_recruitment

def load_qian(figc_path, figd_path):
    # Expect CSV columns at least: Label/Condition, Mean, Error (or SEM)
    dfC = pd.read_csv(figc_path)
    dfD = pd.read_csv(figd_path)
    return dfC, dfD

def extract_targets(dfC, dfD):
    # Fig C: Mets vs Control in Lung (IM:RM)
    # assume rows: Condition: Control Lung, Mets Lung
    def get_mean(df, key):
        row = df[df.iloc[:,0].str.contains(key, case=False, regex=False)]
        return float(row["Mean"].values[0])

    control_lung = get_mean(dfC, "Control")
    mets_lung = get_mean(dfC, "Mets")
    R_lung = mets_lung / control_lung if control_lung>0 else np.nan

    # Fig D: Control Ab vs anti-CCL2 Ab — define blockade ratio (anti / ctrl) < 1
    ctrl_ab = get_mean(dfD, "Ctrl")
    anti_ab = get_mean(dfD, "anti")
    R_block = anti_ab / ctrl_ab if ctrl_ab>0 else np.nan
    return R_lung, R_block

def fit_params(figc_path, figd_path, results_dir):
    dfC, dfD = load_qian(figc_path, figd_path)
    R_lung, R_block = extract_targets(dfC, dfD)

    # Time horizon (Qian often shows ~7h; allow 0..8h)
    times = np.linspace(0, 8.0, 81)

    # Simple loss: difference between model ratios and targets
    def model_ratios(alpha, K_L, dMo, eta, cC, cM2, dL):
        # Control scenario (low tumor drive): set C_drive=0.2, M2_drive=0.0
        pars = dict(alpha=alpha, K_L=K_L, dMo=dMo, cC=cC, cM2=cM2, dL=dL)
        _, Mo_ctrl = simulate_recruitment(times, pars, C_drive=0.2, M2_drive=0.0, eta=1.0, L0=0.1, Mo0=0.1)
        # Mets scenario (tumor present): higher C_drive (+ optional M2 feedback)
        _, Mo_mets = simulate_recruitment(times, pars, C_drive=1.0, M2_drive=0.3, eta=1.0, L0=0.1, Mo0=0.1)
        # Blockade on mets: reduce inflow by eta
        _, Mo_block = simulate_recruitment(times, pars, C_drive=1.0, M2_drive=0.3, eta=eta, L0=0.1, Mo0=0.1)

        # pick final time
        Mc = Mo_ctrl[-1]; Mm = Mo_mets[-1]; Mb = Mo_block[-1]
        Rlung_hat = Mm / Mc if Mc>0 else np.inf
        Rblock_hat = Mb / Mm if Mm>0 else np.inf
        return Rlung_hat, Rblock_hat

    best = None; best_loss = 1e9
    # Crude random/grid search
    rng = np.random.default_rng(42)
    for _ in range(600):
        alpha = 0.2 + 1.8*rng.random()
        K_L  = 0.05 + 1.0*rng.random()
        dMo  = 0.02 + 0.3*rng.random()
        eta  = 0.3 + 0.7*rng.random()  # 0.3..1
        cC   = 0.05 + 0.5*rng.random()
        cM2  = 0.00 + 0.5*rng.random()
        dL   = 0.05 + 0.6*rng.random()
        Rlung_hat, Rblock_hat = model_ratios(alpha,K_L,dMo,eta,cC,cM2,dL)
        loss = (Rlung_hat - R_lung)**2 + (Rblock_hat - R_block)**2
        if loss < best_loss:
            best_loss = loss
            best = dict(alpha=alpha, K_L=K_L, dMo=dMo, eta=eta, cC=cC, cM2=cM2, dL=dL,
                        Rlung_hat=Rlung_hat, Rblock_hat=Rblock_hat)

    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "qian_fit.json"), "w") as f:
        json.dump(best, f, indent=2)

    # Simple bar plot: target vs model
    labels = ["Lung Mets/Control", "anti-CCL2 / Mets"]
    target = [R_lung, R_block]
    model  = [best["Rlung_hat"], best["Rblock_hat"]]
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(6,4))
    ax.bar(x - width/2, target, width, label="Data")
    ax.bar(x + width/2, model,  width, label="Model")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15)
    ax.set_ylabel("Ratio")
    ax.set_title("Qian 2011 — Recruitment ratios")
    ax.legend()
    plt.tight_layout()
    figpath = os.path.join(results_dir, "qian_fit.png")
    plt.savefig(figpath, dpi=160)
    plt.close(fig)
    return best, figpath

if __name__ == "__main__":
    try:
        from config import get_qian_paths
        paths = get_qian_paths()
    except Exception:
        here = os.path.join(os.path.dirname(__file__), "..", "data")
        paths = [os.path.join(here, fn) for fn in ["FigB.csv","FigC.csv","FigD.csv","FigE.csv","FigF.csv"]]
    figc, figd = paths[1], paths[2]  # FigC and FigD as primary for fit
    results_dir = os.path.join(os.path.dirname(__file__), "..", "results", "figures")
    best, figpath = fit_params(figc, figd, results_dir)
    print("Saved:", figpath)
    print("Params:", json.dumps(best, indent=2))
