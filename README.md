# Macrophage–Tumor Interaction Model (new, bone-aware)

**What this does (short):**  
Rebuild of my model after discussing scope with Claire. It keeps things minimal but biologically correct:
- CCL2/CCR2 **recruitment** (Qian),
- TNBC cytokine-driven **drift** to TAM/M2 (Pe),
- myeloid-VEGF **angiogenesis** (Qian),
- **bone awareness** via θ_OC (fraction of monocytes that don’t become macrophages).

**States:**  
C(t): tumor burden (MDA-MB-231)  
A_M2(t): pro-tumor macrophage activity (TAM/M2)  
V(t): angiogenic drive

**Equations:**  
- dC/dt = r·C·(1 − C/K)·(1 + α_V·V) + α_M2·A_M2·C  
- dA_M2/dt = k_drift·S_cyt + k_I·g(C)·σ_CCL2·(1−θ_OC) − d_M2·A_M2  
- dV/dt = s_V·A_M2 − d_V·V  
  with g(C) = C/(C + K_g) and S_cyt = S_scale·(w_IL6·IL6 + w_IL10·IL10 + w_TNF·TNF).

**Data anchors:**  
- Pe Fig 3B → α_M2 from 72h Δr/day (TAM vs control)  
- Pe Fig 2D → IL-6/IL-10/TNF levels for drift  
- Qian Fig 1 → σ_CCL2 (anti-CCL2/control recruitment ratio)  
- Qian Fig 2–5 → κ_VEGF (KO/WT), used to sanity-check angiogenesis

**Run:**  
```bash
python3 -m pip install -r requirements.txt
python3 main.py
