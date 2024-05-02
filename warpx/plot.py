## Energy evolution in the system with your input parameters
import matplotlib.pyplot as plt
from utils.plot import plot_energy_evolution

import json

import os
from pathlib import Path

from utils import export_ts, ts2xr, load_ts_all

import xarray as xr
import numpy as np
import xrft


def normalize_dft_xr(da, meta):
    k_norm = 1 / meta["d_i"]
    w_norm = meta["w_ci"]
    for coord in da.coords:
        if "freq" in coord and "time" in coord:
            da = da.assign_coords({f"{coord}_norm": 2 * np.pi * da[coord] / w_norm})
        if "freq" in coord and "time" not in coord:
            da = da.assign_coords({f"{coord}_norm": 2 * np.pi * da[coord] / k_norm})
    return da


def plot_wk_spectrum_ds(ds: xr.Dataset, fields, meta, step=8):

    k_norm = 1 / meta["d_i"]
    w_norm = meta["w_ci"]

    fig, axes = plt.subplots(len(fields), 2, figsize=(12, 5))

    for i, field in enumerate(fields):
        da = ds[("boxlib", field)]

        # DFT
        da_fft: xr.DataArray = xrft.fft(da, dim=["z", "time"])
        da_fft_mean = da_fft.mean(["x", "y"]).pipe(normalize_dft_xr, meta)
        amp = np.abs(da_fft_mean)
        vmin = -3
        ax0 = axes[i][0]
        ax1 = axes[i][1]

        ax0.pcolormesh(
            amp.freq_z_norm,
            amp.freq_time / w_norm,
            np.log10(amp / amp.max()),
            vmin=vmin,
        )

        # Power spectrum
        p_spec: xr.DataArray = xrft.power_spectrum(da, dim="z").pipe(
            normalize_dft_xr, meta
        )
        ps_mean = p_spec.mean(["x", "y"])
        ps_mean[::step].plot.line(
            x="freq_z_norm",
            hue="time",
            xscale="log",
            yscale="log",
            ylim=(1e-15, 1e-8),
            ax=ax1,
        )

        ax0.set_xlim(-0.4, 0.4)
        ax0.set_xlabel(r"$k d_i$")
        ax0.set_ylabel(r"$\omega / \Omega_i$")
        ax1.set_xlabel(r"$k d_i$")

    return fig, axes


###

from functools import partial
import yt

direction = "z"
ps = "ion"


_field_weight_field = ("boxlib", "volume")
_part_weight_field = (ps, "particle_weight")

plot_field_profile = partial(
    yt.ProfilePlot,
    x_field=direction,
    weight_field=_field_weight_field,
    x_log=False,
    y_log=False,
)

create_field_profile = partial(
    yt.create_profile,
    bin_fields=direction,
    weight_field=_field_weight_field,
    deposition="cic",
)

bin_fields = (ps, f"particle_position_{direction}")

plot_particle_profile = partial(
    yt.ProfilePlot,
    x_field=bin_fields,
    weight_field=_part_weight_field,
    x_log=False,
    y_log=False,
)

create_part_profile = partial(
    yt.create_profile,
    bin_fields=bin_fields,
    weight_field=_part_weight_field,
    deposition="cic",
)


field = "particle_momentum_y"


def plot_field_with_plasma_profile(ds_field, ds_part, field0, field1, twin=True):

    n_bins = ds_field.domain_dimensions[0]

    field_profile = create_field_profile(
        ds_field.all_data(), fields=field0, n_bins=n_bins
    )
    part_profile = create_part_profile(ds_part.all_data(), fields=field1)

    p0 = yt.ProfilePlot.from_profiles(field_profile)
    p1 = yt.ProfilePlot.from_profiles(part_profile)
    p0.set_log(field0, False)
    p1.set_log(field1, False)

    # `ProfilePlot` does not support `deposition`
    # p0 = plot_field(ds_field, y_fields=field0)
    # p1 = plot_particle_profile(ds_part, y_fields=field1)

    # Customizing the plot
    fig, ax = plt.subplots()

    if twin:
        ax2 = ax.twinx()
    else:
        ax2 = ax

    plot = p0.plots[field0]
    plot.figure = fig
    plot.axes = ax

    plot = p1.plots[field1]
    plot.figure = fig
    plot.axes = ax2

    p0.render()
    p1.render()

    if twin:
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.set_offset_position("right")
        # set the color of the labels and lines
        ax2.yaxis.label.set_color("orange")
        ax2.lines[0].set_color("orange")
        # ax2.set_ylim(ax.get_ylim())
        fig.set_size_inches(12, 5)

    ax.set_ylim(-220, 220)
    ax2.set_ylim(-220, 220)
    return fig


def plot_field_with_plasma_profile_ts(
    ts_field,
    ts_part,
    name=None,
    step=8,
    field0=("boxlib", "V_Alfven_y"),
    field1=("V_y"),
):
    os.makedirs("figures/field_plasma_profile", exist_ok=True)
    for ds_field, ds_part in zip(ts_field[::step], ts_part[::step]):
        add_field(ds_field)
        add_field(ds_part)
        fig = plot_field_with_plasma_profile(ds_field, ds_part, field0, field1)
        fig.savefig(f"figures/field_plasma_profile/{name}_{ds_field.current_time}.png")
        plt.close(fig)


###

import typer

app = typer.Typer()


@app.command()
def main(
    dim: int = 1,
    beta: float = 0.25,
    theta: float = 60,
    plasma_resistivity: float = 100,
    export: bool = False,
    plot_wk_spectrum: bool = False,
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


def _plot_field_with_plasma_profile(ds_field, ds_part, field0, field1, twin=False):
    # BUG: not working, from_profiles does not support different data sources

    field_profile = create_field_profile(ds_field.all_data(), fields=field0)
    part_profile = create_part_profile(ds_part.all_data(), fields=field1)

    field_label = field0[-1]
    part_label = field1[-1]

    p = yt.ProfilePlot.from_profiles(
        [field_profile, part_profile],
        labels=[field_label, part_label],
    )

    return p.figure


if __name__ == "__main__":
    app()
