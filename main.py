
# main.py — end-to-end: run Pe & Qian fits, then a coupled demo plot
import os, json
import numpy as np
import matplotlib.pyplot as plt

from fit.fit_pe import fit_params as fit_pe
from fit.fit_qian import fit_params as fit_qian
from model.equations import simulate_recruitment, simulate_proliferation

BASE = os.path.dirname(__file__)
RES_DIR = os.path.join(BASE, "results", "figures")
os.makedirs(RES_DIR, exist_ok=True)

# 1) Fit Pe proliferation
try:
    from config import get_pe_path
    pe_path = get_pe_path()
except Exception:
    pe_path = os.path.join(BASE, "data", "Pe_Proliferation.csv")
pe_params, pe_fig = fit_pe(pe_path, RES_DIR)

# 2) Fit Qian recruitment
try:
    from config import get_qian_paths
    paths = get_qian_paths()
except Exception:
    here = os.path.join(BASE, "data")
    paths = [os.path.join(here, fn) for fn in ["FigB.csv","FigC.csv","FigD.csv","FigE.csv","FigF.csv"]]
figc, figd = paths[1], paths[2]
qian_best, qian_fig = fit_qian(figc, figd, RES_DIR)

# Save combined params
combined = {"pe": pe_params, "qian": qian_best}
with open(os.path.join(BASE, "results", "combined_params.json"), "w") as f:
    json.dump(combined, f, indent=2)

# 3) Coupled demo: show C(t) growth under scenarios using pe params
times = np.array([0, 24, 48, 72], dtype=float)
parsP = dict(r0=pe_params["r0"], gamma2=pe_params["gamma2"], K2=pe_params["K2"],
             gamma1=pe_params["gamma1"], K1=pe_params["K1"])

ctrl = dict(S_M1=0.0, S_M2=0.0, C0=1.0)
M2   = dict(S_M1=0.0, S_M2=pe_params["s_M2"], C0=1.0)
TAM  = dict(S_M1=pe_params["theta"]*pe_params["s_TAM"],
            S_M2=(1.0-pe_params["theta"])*pe_params["s_TAM"], C0=1.0)

C_ctrl = simulate_proliferation(times, parsP, ctrl)
C_m2   = simulate_proliferation(times, parsP, M2)
C_tam  = simulate_proliferation(times, parsP, TAM)

fig, ax = plt.subplots(figsize=(6,4))
ax.plot(times, np.ones_like(times), label="Control (model)")           # << fixed: always 1
ax.plot(times, C_m2 / C_ctrl,      label="+M2 (model)")
ax.plot(times, C_tam / C_ctrl,     label="TAM (model)")
ax.set_xlabel("Time (h)")
ax.set_ylabel("Relative fold vs Control")
ax.set_title("Coupled demo — Proliferation scenarios")
ax.legend()
plt.tight_layout()
demo_fig = os.path.join(RES_DIR, "coupled_demo.png")
plt.savefig(demo_fig, dpi=160)
plt.close(fig)


print("Saved figures:")
print("-", pe_fig)
print("-", qian_fig)
print("-", demo_fig)
