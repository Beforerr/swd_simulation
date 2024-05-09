import os
from pathlib import Path

def change_dir(dim, beta, theta, eta, wave_length, base_dir: Path = Path().cwd()):
    sub_dir = f"dim_{dim}_beta_{beta}_theta_{theta}_eta_{eta}_l_{wave_length}"
    directory = base_dir / sub_dir
    os.makedirs(directory, exist_ok=True)
    os.chdir(directory)
