import matplotlib.pyplot as plt
import numpy as np

def plot_timecourses(ts, ys, title=None):
    C = ys[:, 0]
    A = ys[:, 1]
    V = ys[:, 2]

    fig, ax1 = plt.subplots(figsize=(7, 4))

    # Left axis: fold-change of C (so we see growth shape clearly)
    C_fold = C / max(C[0], 1e-12)
    ax1.plot(ts, C_fold, label="C (fold)", linewidth=2)
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Tumor C / C0 (fold)")
    if title: ax1.set_title(title)

    # Right axis: A_M2 and V on their natural scales
    ax2 = ax1.twinx()
    ax2.plot(ts, A, label="A_M2", linestyle="--")
    ax2.plot(ts, V, label="V", linestyle=":")
    ax2.set_ylabel("A_M2, V (a.u.)")

    # Build a joint legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    fig.tight_layout()
    return fig, (ax1, ax2)

def plot_scenarios(ts_list, ys_list, names, title="C(t) scenarios"):
    fig, ax = plt.subplots(figsize=(7, 4))
    for ts, ys, nm in zip(ts_list, ys_list, names):
        C_fold = ys[:, 0] / max(ys[0, 0], 1e-12)
        ax.plot(ts, C_fold, label=nm, linewidth=2)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Tumor C / C0 (fold)")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    return fig, ax
