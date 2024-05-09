import os
from pathlib import Path
import json

def setup(dim, beta, theta, eta, base_dir = Path("oblique_linear_alfven")):
    # change to simulation directory and load metadata
    sub_dir = f"dim_{dim}_beta_{beta}_theta_{theta}_eta_{eta}"
    directory = base_dir / sub_dir

    try:
        os.makedirs(directory, exist_ok=True)
        os.chdir(directory)
        os.makedirs("figures", exist_ok=True)
    finally:
        pass

    # load simulation metadata (json)
    with open("sim_parameters.json") as f:
        data = json.load(f)
    return data