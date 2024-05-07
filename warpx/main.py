## Energy evolution in the system with your input parameters
import matplotlib.pyplot as plt
import json
from utils import export_ts, ts2xr, load_ts_all
from utils.plot import plot_wk_spectrum_ds
from utils.energy import plot_energy_evolution
from utils.pressure import export_pressure_field

import os
from pathlib import Path

import typer

app = typer.Typer()


@app.command()
def main(
    dim: int = 1,
    beta: float = 0.25,
    theta: float = 60,
    plasma_resistivity: float = 100,
    export: bool = True,
    export_pressure: bool = False,
    plot_wk_spectrum: bool = True,
):

    base_dir = Path(os.getcwd()) / "01_oblique_linear_alfven"
    sub_dir = f"dim_{dim}_beta_{beta}_theta_{theta}_eta_{plasma_resistivity}"
    directory = base_dir / sub_dir
    os.chdir(directory)
    os.makedirs("figures", exist_ok=True)

    # load simulation parameters
    with open("sim_parameters.json", "rb") as f:
        meta = json.load(f)

    plot_energy_evolution(meta=meta)
    plt.savefig("figures/energy_evolution.png")

    ts_field, ts_part = load_ts_all(meta)

    if export:
        print("Exporting data...")
        export_ts(ts_field)

    if plot_wk_spectrum:
        print("Plotting wk spectrum...")
        fields = ["By", "Bx"]
        ds = ts2xr(ts_field)
        plot_wk_spectrum_ds(ds, meta=meta, fields=fields)
        plt.savefig("figures/wk_spectrum.png")
        plt.savefig("figures/wk_spectrum.pdf")
        
    if export or export_pressure:
        print("Exporting pressure field...")
        export_pressure_field(ts_part, ts_field, meta=meta)

if __name__ == "__main__":
    app()
