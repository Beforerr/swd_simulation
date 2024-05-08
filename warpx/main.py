## Energy evolution in the system with your input parameters
import matplotlib.pyplot as plt
import json
from utils import export_ts, ts2xr, load_ts_all
from utils.plot import (
    plot_wk_spectrum_ds,
    plot_field_with_plasma_profile_ts,
    plot_energy_evolution,
)
from utils.pressure import export_pressure_field

import os
from pathlib import Path

import typer

app = typer.Typer()

def setup(dim, beta, theta, eta):

    # change to simulation directory and load metadata
    base_dir = Path(os.getcwd()) / "01_oblique_linear_alfven"
    sub_dir = f"dim_{dim}_beta_{beta}_theta_{theta}_eta_{eta}"
    directory = base_dir / sub_dir

    try:
        os.chdir(directory)
        os.makedirs("figures", exist_ok=True)
    except:
        pass

    # load simulation metadata (json)
    with open("sim_parameters.json") as f:
        data = json.load(f)
    return data


@app.command()
def main(
    dim: int = 1,
    beta: float = 0.25,
    theta: float = 60,
    eta: float = 100,
    export: bool = True,
    plot_energy = True,
    plot_wk_spectrum: bool = False,
    plot_va_vl: bool = False
):

    meta = setup(dim, beta, theta, eta)

    ts_field, ts_part = load_ts_all(meta)

    if export:
        print("Exporting data...")
        export_ts(ts_field)
        print("Exporting pressure field...")
        export_pressure_field(ts_part, ts_field, meta=meta)

    if plot_energy:
        plot_energy_evolution(meta=meta)
        plt.savefig("figures/energy_evolution.png")

    if plot_wk_spectrum:
        print("Plotting wk spectrum...")
        fields = ["By", "Bx"]
        ds = ts2xr(ts_field)
        plot_wk_spectrum_ds(ds, meta=meta, fields=fields)
        plt.savefig("figures/wk_spectrum.png")
        plt.savefig("figures/wk_spectrum.pdf")
        
    if plot_va_vl:
        plot_field_with_plasma_profile_ts(
            ts_field=ts_field, ts_part=ts_part,meta=meta,
            name="va_vl"
        )


if __name__ == "__main__":
    app()
