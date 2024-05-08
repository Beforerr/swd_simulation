## Energy evolution in the system with your input parameters
import matplotlib.pyplot as plt
import os

import xarray as xr
import numpy as np
import xrft

from .energy import plot_energy_evolution
from .yt_warpx import add_field
from matplotlib.axes import Axes
from tqdm import tqdm


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


###

from functools import partial
import yt

direction = "z"
ps = "ions"

_field_weight_field = ("boxlib", "volume")
_part_weight_field = (ps, "particle_weight")

plot_field_profile = partial(
    yt.ProfilePlot,
    x_field=direction,
    weight_field=_field_weight_field,
    x_log=False,
    y_log=False,
)


def plot_field_with_plasma_profile(ds_field, ds_part, field0, field1, twin=True):

    n_bins = ds_field.domain_dimensions[0]

    if ds_part.dimensionality == 1:
        direction = "x"
        true_direction = "z"

    field_profile = yt.create_profile(
        ds_field.all_data(),
        fields=field0,
        n_bins=n_bins,
        bin_fields=direction,
        weight_field=_field_weight_field,
        deposition="cic",
    )

    _part_bin_field = (ps, f"particle_position_{direction}")
    part_profile = yt.create_profile(
        ds_part.all_data(),
        fields=field1,
        bin_fields=_part_bin_field,
        weight_field=_part_weight_field,
        deposition="cic",
    )

    # NOTE: from_profiles does not support different data sources
    p0 = yt.ProfilePlot.from_profiles(field_profile)
    p1 = yt.ProfilePlot.from_profiles(part_profile)
    p0.set_log(field0, False)
    p1.set_log(field1, False)

    # `ProfilePlot` does not support `deposition`
    # p0 = plot_field(ds_field, y_fields=field0)
    # p1 = plot_particle_profile(ds_part, y_fields=field1)

    # Customizing the plot
    fig, ax = plt.subplots()
    ax: Axes

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


def plot_plasma_velocity_profile_ts(
    ts_part,
    step=8,
    x_bins=128,
    y_bins=64
):
    directory = "figures/plasma_velocity/"
    os.makedirs(directory, exist_ok=True)
    xfield = ("ions", "particle_position_x")
    yfields = [
        ("ions", "particle_momentum_y"),
        ("ions", "particle_momentum_z"),
        ("ions", "particle_momentum_x"),
    ]

    def plot(yfield):
        p = yt.ParticlePlot(
            ds_part,
            xfield,
            yfield,
            z_fields=_part_weight_field,
            x_bins=x_bins,
            y_bins=y_bins
        )
        p.set_log(field="all", log=False)
        p.save(directory)
        return p

    for ds_part in tqdm(ts_part[::step]):
        list(map(plot, yfields))


def plot_field_with_plasma_profile_ts(
    ts_field,
    ts_part,
    meta,
    name=None,
    step=8,
    field0=("boxlib", "V_Alfven_y"),
    field1=("V_y"),
):
    directory = "figures/field_plasma_profile"
    os.makedirs(directory, exist_ok=True)
    for ds_field, ds_part in zip(ts_field[::step], ts_part[::step]):
        add_field(ds_field, meta=meta)
        add_field(ds_part, meta=meta)
        fig = plot_field_with_plasma_profile(ds_field, ds_part, field0, field1)

        fname = f"{name}_{round(ds_field.current_time.item())}.png"
        fig.savefig(f"{directory}/{fname}")
        plt.close(fig)


def hodogram_ds(ds, meta: dict, comp1="By", comp2="Bz"):
    time = ds.current_time
    time_norm = time.value / meta["t_ci"]
    ad = ds.all_data()
    plt.plot(ad[comp1], ad[comp2], label=f"t={time_norm:.2f}")
    plt.xlabel(comp1)
    plt.ylabel(comp2)
