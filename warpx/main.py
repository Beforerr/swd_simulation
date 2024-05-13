## Energy evolution in the system with your input parameters
import matplotlib.pyplot as plt
from utils import setup
from utils.io import export_ts, ts2xr, load_ts_all
from utils.plot import (
    plot_wk_spectrum_ds,
    plot_field_with_plasma_profile_ts,
    plot_energy_evolution,
    plot_plasma_velocity_profile_ts,
)
from utils.pressure import export_part_field
import typer

def info(meta):
    from plasmapy.formulary import gyroradius
    import astropy.units as u
    B = meta["B0"] * u.T
    T = meta["T_plasma"] * u.eV
    
    r = gyroradius(B, particle="p+", T=T)
    print(f"Gyroradius: {r.to(u.km)}")
    
app = typer.Typer()

@app.command()
def main(
    directory: str = ".",
    export: bool = True,
    plot_energy: bool = True,
    plot_wk_spectrum: bool = False,
    plot_va_vl: bool = False,
    plot_plasma_velocity: bool = False
):

    meta = setup(directory)
    info(meta)

    ts_field, ts_part = load_ts_all(meta)

    if export:
        print("Exporting data...")
        export_ts(ts_field)
        print("Exporting pressure field...")
        export_part_field(ts_part, ts_field, meta=meta)

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

    if plot_plasma_velocity:
        print("Plotting plasma velocity profile...")
        plot_plasma_velocity_profile_ts(ts_part=ts_part)

if __name__ == "__main__":
    app()
