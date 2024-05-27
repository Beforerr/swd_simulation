import os
import json
from beforerr.project import datadir


def setup(sub_dir, base_dir=datadir()):
    # change to simulation directory and load metadata
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
