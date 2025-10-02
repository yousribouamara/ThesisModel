def pe_condition(condition: str, theta=0.3):
    # Map Pe conditions to scenario dict for the proliferation fit.
    # Control: no macrophage stimulus
    # M2:     boosted M2 activity
    # TAM:    mixed M1/M2 (theta controls the mix)
    scen = dict(use_ratio_shortcut=True, theta=theta, blockade_eta=1.0)
    condition_l = condition.lower()
    if condition_l == "control":
        scen.update(C0=1.0, L0=0.1, Mo0=0.1)
    elif condition_l == "m2":
        scen.update(C0=1.0, L0=0.2, Mo0=0.15)
    elif condition_l == "tam":
        scen.update(C0=1.0, L0=0.2, Mo0=0.15, theta=theta)
    else:
        raise ValueError(f"Unknown Pe condition: {condition}")
    return scen

def qian_scenario(kind: str, blockade=False):
    # Build Qian-like scenarios:
    # - kind: 'control_lung' vs 'mets_lung'
    # - blockade: True applies eta < 1
    eta = 0.5 if blockade else 1.0
    if kind == "control_lung":
        return dict(use_ratio_shortcut=True, theta=0.6, L0=0.1, C0=0.5, Mo0=0.1, blockade_eta=eta)
    elif kind == "mets_lung":
        # higher tumor drive → higher L
        return dict(use_ratio_shortcut=True, theta=0.4, L0=0.3, C0=1.0, Mo0=0.1, blockade_eta=eta)
    else:
        raise ValueError(f"Unknown Qian scenario: {kind}")
