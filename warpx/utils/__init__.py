import os
import json

def setup(directory):
    # change to simulation directory and load metadata
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