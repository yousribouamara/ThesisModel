from dataclasses import dataclass
from typing import Dict, Tuple
import numpy as np

@dataclass
class Params:
    r: float = 0.8
    K: float = 1e9
    alpha_M2: float = 0.2
    k_drift: float = 0.8
    d_M2: float = 0.3
    w_IL6: float = 1.0
    w_IL10: float = 0.5
    w_TNF: float = 0.2
    S_scale: float = 1e-3
    s_V: float = 0.05
    d_V: float = 0.2
    alpha_V: float = 0.5
    vegf_enabled: bool = True
    k_I: float = 1.0
    K_g: float = 1e6
    sigma_CCL2: float = 1.0
    theta_OC: float = 0.2

def _g_size(C, K_g): return C / (C + K_g) if (C + K_g) > 0 else 0.0

def rhs(t, y, p: Params, drivers: Dict[str, float]):
    C, A_M2, V = y
    I = p.k_I * _g_size(C, p.K_g) * p.sigma_CCL2 * (1 - p.theta_OC)
    S = p.S_scale * (p.w_IL6*drivers.get("IL6", 0) + p.w_IL10*drivers.get("IL10", 0) + p.w_TNF*drivers.get("TNF", 0))
    dA_M2 = p.k_drift * S + I - p.d_M2 * A_M2
    sV = p.s_V if p.vegf_enabled else 0
    dV = sV * A_M2 - p.d_V * V
    growth = p.r * C * (1 - C/p.K) * (1 + p.alpha_V * V)
    dC = growth + p.alpha_M2 * A_M2 * C
    return np.array([dC, dA_M2, dV])

def simulate(tspan: Tuple[float,float], y0: Tuple[float,float,float],
             p: Params, drivers: Dict[str,float], dt=0.01):
    t0, t1 = tspan
    n = int(np.ceil((t1 - t0) / dt)) + 1
    ts = np.linspace(t0, t1, n)
    ys = np.zeros((n, 3))
    ys[0] = y0
    for i in range(n - 1):
        t, y = ts[i], ys[i]
        k1 = rhs(t, y, p, drivers)
        k2 = rhs(t + dt/2, y + dt*k1/2, p, drivers)
        k3 = rhs(t + dt/2, y + dt*k2/2, p, drivers)
        k4 = rhs(t + dt, y + dt*k3, p, drivers)
        ys[i+1] = np.maximum(y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4), 0)
    return ts, ys
