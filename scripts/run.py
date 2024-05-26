import typer
from rich import print
from beforerr.project import setup_run_dir
from warpx.oblique_linear_alfven import AlfvenModes
app = typer.Typer()

@app.command(
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True}
)
def main(
    ctx: typer.Context,
    dim: int = 1,
    beta: float = 0.25,
    theta: float = 60,
    eta: float = 100,
    wave_length: float = 64,
    Te_norm: float = 1,
    dz_norm: float = 0.5,
    dt_norm: float = 1 / 64 ,
    Lz_norm: float = 64,
    time_norm: float = 100,
    substeps: int = 16,
    nppc: int = 64,
    dry_run: bool = False,
    verbose: bool = False,
):

    setup_run_dir(ctx.params, accesses=["dim", "beta", "theta", "eta", "wave_length", "Te_norm"])

    sim_kwargs = dict()
    if dim == 3:
        # Reduce the number of cells in x and y to accelerate the simulation
        sim_kwargs.update(
            nx=8,
            ny=8,
            grid_kwargs=dict(
                warpx_blocking_factor_x=4,
                warpx_blocking_factor_y=4,
            ),
        )

    simulation = AlfvenModes(
        **ctx.params,
        plasma_resistivity=eta,
        diag_part=True,
        **sim_kwargs
    )

    print(ctx.params)
    # wavenumber
    print("wavenumer (1 / ion inertial length):", simulation.k * simulation.d_i)

    # wave period
    w = simulation.k * simulation.vA
    t_w = 2 * np.pi / w
    t_w_norm = t_w / simulation.t_ci
    t_w_norm

    simulation._sim.verbose = verbose
    dry_run or simulation._sim.step()
    
    return simulation


if __name__ == "__main__":
    app()