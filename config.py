

"""
# Macrophage–Tumor Interaction Model (Fresh Build)

Author: Yousri Bouamara
Minimal, clean, bone-aware model based on Pe et al. and Qian et al.

**Equations**

- dC/dt = r·C·(1 − C/K)·(1 + α_V·V) + α_M2·A_M2·C
- dA_M2/dt = k_drift·S_cyt + k_I·g(C)·σ_CCL2·(1 − θ_OC) − d_M2·A_M2
- dV/dt = s_V·A_M2 − d_V·V


Note for Claire:
Central place for file paths. I pointed this at my local Mac folders.
If I move the data, I’ll just update these strings.
"""

from pathlib import Path

# Where my digitized CSVs are (Pe & Qian)
# For Claire: change these to where you save them '.../Downloads/DataPe'. I'll add a zip in mail
DATA_PE   = Path("/Users/yousribouamara/Documents/Thesis/Data/DataPe")
DATA_QIAN = Path("/Users/yousribouamara/Documents/Thesis/Data/DataQian")

# Output folder for figures / summaries
RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True, parents=True)






