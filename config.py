
# config.py — centralizes data paths (local machine first, then bundled /data fallback)
import os

# User-local default paths (edit these to match your machine)
USER_DATA_DIR_PE = "/Users/yousribouamara/Downloads/ThesisData/DataPe"
USER_DATA_DIR_QIAN = "/Users/yousribouamara/Downloads/ThesisData/DataQian"

# Bundled fallback (inside the project)
BUNDLED_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

def path_or_fallback(local_path, fallback_filename):
    if local_path and os.path.exists(local_path):
        return local_path
    fb = os.path.join(BUNDLED_DATA_DIR, fallback_filename)
    if os.path.exists(fb):
        return fb
    raise FileNotFoundError(f"Neither local path '{local_path}' nor bundled data '{fallback_filename}' found.")

def get_pe_path():
    local = os.path.join(USER_DATA_DIR_PE, "Pe_Proliferation.csv")
    return path_or_fallback(local, "Pe_Proliferation.csv")

def get_qian_paths():
    files = ["FigB.csv","FigC.csv","FigD.csv","FigE.csv","FigF.csv"]
    out = []
    for fn in files:
        local = os.path.join(USER_DATA_DIR_QIAN, fn)
        out.append(path_or_fallback(local, fn))
    return out
