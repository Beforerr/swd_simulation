from rich import print
from beforerr.project import setup_run_dir
from warpx.oblique_linear_alfven import AlfvenModes

from sacred import Experiment

ex = Experiment()

DefaultParams = dict(
    dim=1, beta=0.25, theta=60, plasma_resistivity=100, wave_length=64, Te_norm=1
)
DefaultOtherParams = dict(
    dz_norm=0.5, dt_norm=1 / 64, Lz_norm=64, time_norm=100, substeps=16, nppc=64
)

ex.add_config(
    params=DefaultParams,
    params_s=DefaultOtherParams,
)


def setup_sim_params(params):
    sim_kwargs = dict()
    grid_kwargs = dict(warpx_blocking_factor_x=4, warpx_blocking_factor_y=4)
    if params.get("dim") == 3:
        # Reduce the number of cells in x and y to accelerate the simulation
        sim_kwargs.update(nx=8, ny=8, grid_kwargs=grid_kwargs)
    return sim_kwargs


@ex.automain
def main(
    params: dict,
    params_s: dict,
    dry_run: bool = False,
    verbose: bool = False,
):
    setup_run_dir(params, sort=False)
    sim_kwargs = setup_sim_params(params)

    simulation = AlfvenModes(**params, **params_s, **sim_kwargs, diag_part=True)

    # wavenumber
    print("wavenumer (1 / ion inertial length):", simulation.k * simulation.d_i)

    simulation._sim.verbose = verbose
    dry_run or simulation._sim.step()

    return simulation