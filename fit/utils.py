
# fit/utils.py — small helpers for weighted least squares without scipy
import numpy as np

def wls_loss(y_pred, y_obs, w=None):
    y_pred = np.asarray(y_pred, float)
    y_obs = np.asarray(y_obs, float)
    if w is None:
        w = np.ones_like(y_obs)
    else:
        w = np.asarray(w, float)
    return np.sum(w * (y_pred - y_obs)**2)

def sem_to_weight(sem, eps=1e-9):
    sem = np.asarray(sem, float)
    return 1.0 / (sem**2 + eps)

def bounded(val, lo, hi):
    return max(lo, min(hi, val))

def grid_search_1d(func, lo, hi, n=50):
    xs = np.linspace(lo, hi, n)
    vals = [func(x) for x in xs]
    i = int(np.argmin(vals))
    return xs[i], vals[i]
