
## Energy evolution in the system with your input parameters
import matplotlib.pyplot as plt
from utils.plot import plot_energy_evolution

plot_energy_evolution(meta=meta)
plt.savefig("figures/energy_evolution.png")

from utils import ds2df, export_ts, ts2xr

export_ts(ts_field)

ds = ts2xr(ts_field)
ds.assign_coords(time_norm = ds.time / meta["t_ci"])

import xrft
import xarray as xr

def plot_wk_spectrum(ds: xr.Dataset, fields, step=8):
    fig, axes = plt.subplots(len(fields), 2, figsize=(12, 5))

    for i, field in enumerate(fields):
        da = ds[("boxlib", field)]

        # DFT
        da_fft: xr.DataArray = xrft.fft(da, dim=["z", "time"])
        da_fft_mean = da_fft.mean(["x", "y"]).pipe(normalize_dft_xr)
        amp = np.abs(da_fft_mean)

        vmin = -3
        ax0 = axes[i][0]
        ax1 = axes[i][1]
        
        ax0.pcolormesh(
            amp.freq_z_norm, amp.freq_time / w_norm, np.log10(amp / amp.max()), vmin=vmin
        )

        # Power spectrum
        p_spec: xr.DataArray = xrft.power_spectrum(da, dim="z").pipe(normalize_dft_xr)
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

fields = ["By", "Bx"]        
plot_wk_spectrum(ds, fields, step = 16)