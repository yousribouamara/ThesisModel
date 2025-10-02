import numpy as np
from scipy.integrate import solve_ivp
from .equations import rhs_pure, initial_state

def simulate(params, scen, t_end, n_points=200):
    t_span = (0.0, float(t_end))
    y0 = initial_state(params, scen)
    sol = solve_ivp(lambda t, y: rhs_pure(t, y, params, scen),
                    t_span, y0, t_eval=np.linspace(*t_span, n_points),
                    method="RK45", rtol=1e-6, atol=1e-9)
    return sol.t, sol.y  # t, [C, Mo, M1, M2, L]
