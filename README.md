
# ThesisModel

A minimal, modular ODE model of breast-cancer–macrophage interactions calibrated to:
- **Pe et al. 2022** (PLOS ONE): MDA-MB-231 proliferation with TAM/M2 (WST at 24/48/72 h)
- **Qian et al. 2011** (Nature): CCL2-driven monocyte recruitment and anti-CCL2 blockade

## Quick start

```bash
# (optional) python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python main.py
```

Figures land in `results/figures/`:
- `pe_fit.png` — model vs Pe relative fold-changes
- `qian_fit.png` — model vs Qian recruitment ratios
- `coupled_demo.png` — illustrative proliferation scenarios

## Local data paths

Edit `config.py` if your data live elsewhere. By default, it tries your paths:

- /Users/yousribouamara/Downloads/ThesisData/DataPe/Pe_Proliferation.csv
- /Users/yousribouamara/Downloads/ThesisData/DataQian/Fig*.csv

and falls back to the bundled copies in `data/`.

## Structure

- `model/equations.py` — ODE right-hand sides and simple integrators
- `fit/fit_pe.py` — fits r0, gamma2, K2 (and scales) to Pe WST data
- `fit/fit_qian.py` — fits alpha, K_L, dMo, eta (and cC,cM2,dL coarse) to Qian ratios
- `main.py` — runs both fits and produces an illustrative coupled plot
