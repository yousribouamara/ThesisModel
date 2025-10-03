
# model/equations.py — defines ODE right-hand-sides and helpers
# overview on ppt
import numpy as np

def rhs_proliferation(t, C, params, S_M1=0.0, S_M2=0.0):
    """
    dC/dt = r0*C * [ 1 + gamma2 * S_M2/(S_M2+K2) - gamma1 * S_M1/(S_M1+K1) ]
    """
    r0 = params["r0"]
    gamma2 = params.get("gamma2", 0.0)
    gamma1 = params.get("gamma1", 0.0)
    K2 = params.get("K2", 1.0)
    K1 = params.get("K1", 1.0)
    stim = 1.0 + gamma2 * (S_M2 / (S_M2 + K2 + 1e-12)) - gamma1 * (S_M1 / (S_M1 + K1 + 1e-12))
    return r0 * C * stim

def rhs_recruitment(t, Mo, L, params, eta=1.0):
    """
    dMo/dt = alpha * L/(K_L + L) * eta - dMo * Mo
    """
    alpha = params["alpha"]
    K_L = params["K_L"]
    dMo = params.get("dMo", 0.1)
    inflow = alpha * (L / (K_L + L + 1e-12)) * eta
    return inflow - dMo * Mo

def rhs_ccl2(t, L, C, M2, params):
    """
    dL/dt = cC*C + cM2*M2 - dL*L
    """
    cC = params.get("cC", 0.1)
    cM2 = params.get("cM2", 0.0)
    dL = params.get("dL", 0.2)
    return cC * C + cM2 * M2 - dL * L

def simulate_proliferation(times_h, params, scenario):
    """
    Integrate proliferation ODE with explicit Euler (sufficient for short spans).
    scenario: dict with keys S_M1, S_M2, C0
    """
    dt = np.diff(times_h)
    C = np.zeros_like(times_h, dtype=float)
    C[0] = scenario.get("C0", 1.0)
    S_M1 = scenario.get("S_M1", 0.0)
    S_M2 = scenario.get("S_M2", 0.0)
    for i in range(1, len(times_h)):
        t = times_h[i-1]
        dC = rhs_proliferation(t, C[i-1], params, S_M1, S_M2)
        C[i] = max(1e-12, C[i-1] + dt[i-1] * dC)
    return C

def simulate_recruitment(times_h, params, C_drive=1.0, M2_drive=0.0, eta=1.0, L0=0.1, Mo0=0.1):
    """
    Simple coupled Euler integration for L(t) and Mo(t) on a short horizon.
    """
    dt = np.diff(times_h)
    L = np.zeros_like(times_h, dtype=float); L[0] = L0
    Mo = np.zeros_like(times_h, dtype=float); Mo[0] = Mo0
    for i in range(1, len(times_h)):
        t = times_h[i-1]
        # treat C and M2 as constants drives on this short window
        dL = rhs_ccl2(t, L[i-1], C_drive, M2_drive, params)
        L[i] = max(0.0, L[i-1] + dt[i-1]*dL)
        dMo = rhs_recruitment(t, Mo[i-1], L[i], params, eta=eta)
        Mo[i] = max(0.0, Mo[i-1] + dt[i-1]*dMo)
    return L, Mo
