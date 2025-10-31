import numpy as np
import pandas as pd
from config import DATA_PE, DATA_QIAN, RESULTS
from model.equations import Params, simulate
from model.ingest import discover_csvs, parse_pe_3b, parse_pe_2d, parse_qian_sigma_ccl2, parse_qian_kappa_vegf
from model.plotting import plot_timecourses, plot_scenarios

def main():
    # --- discover CSVs ---
    pe_map = discover_csvs(DATA_PE)
    qian_map = discover_csvs(DATA_QIAN)
    print("Pe CSVs:", list(pe_map.keys()))
    print("Qian CSVs:", list(qian_map.keys()))

    # --- extract numeric values from literature ---
    pe3b = parse_pe_3b(pe_map.get("3BPe", DATA_PE / "3BPe.csv"))
    pe2d = parse_pe_2d(pe_map.get("2DPe", DATA_PE / "2DPe.csv"))

    sigma_vals = [parse_qian_sigma_ccl2(qian_map.get(k, DATA_QIAN / f"{k}.csv"))
                  for k in ["Fig1BQian","Fig1CQian","Fig1DQian","Fig1EQian","Fig1FQian"]]
    sigma_vals = [s for s in sigma_vals if s is not None]
    sigma_ccl2 = float(np.mean(sigma_vals)) if sigma_vals else 1.0

    kappa_vegf = parse_qian_kappa_vegf(qian_map.get("Fig2-5Qian", DATA_QIAN / "Fig2-5Qian.csv"))

    drivers_tam = {
        "IL6": pe2d["TAM"].get("IL6", 0.0),
        "IL10": pe2d["TAM"].get("IL10", 0.0),
        "TNF": pe2d["TAM"].get("TNF", 0.0),
    }

    # --- parameter sets ---
    # --- parameter sets ---
    p_ctrl = Params(
        # if you have a 72h control fold from Pe_Proliferation.csv, set r = ln(f_ctrl)/3
        r=0.7,  # slightly smaller baseline growth than 0.8
        K=1e8,  # smaller K so logistic term matters a bit
        alpha_M2=pe3b.get("delta_r_TAM_per_day", 0.25),  # modestly stronger TAM push

        # make recruitment visibly drive A_M2
        k_I=2.5,  # ↑ influx strength
        K_g=2e5,  # ↓ saturation → more sensitive to C
        sigma_CCL2=1.0,
        theta_OC=0.2,

        # drift weights (keep as is; recruitment already adds to A_M2)
        k_drift=0.8, d_M2=0.3, w_IL6=1.0, w_IL10=0.5, w_TNF=0.2, S_scale=1e-3,

        # angiogenesis visible but not crazy
        s_V=0.15,  # ↑ production of V from A_M2
        d_V=0.25,  # slightly faster decay
        alpha_V=0.9,  # V amplifies growth noticeably
        vegf_enabled=True
    )

    # Variants
    p_ccl2 = Params(**p_ctrl.__dict__)
    p_ccl2.sigma_CCL2 = sigma_ccl2  # from Qian; typically < 1

    p_vegf = Params(**p_ctrl.__dict__)
    p_vegf.vegf_enabled = False

    p_bone = Params(**p_ctrl.__dict__)
    p_bone.theta_OC = 0.5  # stronger diversion to osteoclast lineage in bone


    #sanity check
    print("\n[Sanity] Parameters used (control):")
    for k, v in p_ctrl.__dict__.items():
        if isinstance(v, (int, float, bool)):
            print(f"  {k}: {v}")

    print("\n[Sanity] Data-derived:")
    print(f"  alpha_M2 (Pe 3B Δr/day): {p_ctrl.alpha_M2}")
    print(f"  sigma_CCL2 (Qian mean):  {sigma_ccl2}")
    print(f"  VEGF enabled:            {p_ctrl.vegf_enabled}")

    # --- simulation ---
    C0, A20, V0 = 1e5, 0.0, 0.0
    tspan = (0.0, 3.0)  # 72 hours
    runs = []
    for name, P in [("Control", p_ctrl), ("anti-CCL2", p_ccl2), ("VEGF-KO", p_vegf), ("Bone-aware", p_bone)]:
        ts, ys = simulate(tspan, (C0, A20, V0), P, drivers_tam)
        runs.append((name, ts, ys))

    ts_list = [ts for _, ts, _ in runs]
    ys_list = [ys for _, _, ys in runs]
    names   = [nm for nm, _, _ in runs]

    fig1, _ = plot_scenarios(ts_list, ys_list, names)
    fig1.savefig(RESULTS / "Ct_scenarios.png", dpi=200)
    fig2, _ = plot_timecourses(runs[0][1], runs[0][2])
    fig2.savefig(RESULTS / "control_timecourses.png", dpi=200)

    pd.DataFrame([pe3b]).to_csv(RESULTS / "Pe3B_summary.csv", index=False)
    pd.DataFrame([{"sigma_CCL2_mean": sigma_ccl2, "kappa_VEGF": kappa_vegf}]).to_csv(RESULTS / "Qian_summary.csv", index=False)

    print("Results saved in 'results/'")

if __name__ == "__main__":
    main()
