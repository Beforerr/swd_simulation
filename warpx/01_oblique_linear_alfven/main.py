import numpy as np
from pywarpx import picmi
from space_analysis.utils.math import cosd, sind
from space_analysis.simulation.warpx import HybridSimulation
import os
from pathlib import Path
import typer
from rich import print

constants = picmi.constants

def init_obl_alfven(
    k,
    B0,
    A,  #: relative amplitude
    theta,
    vA,
    n0,
    v_ti,  #: ion thermal velocity
):
    """
    Generate field and distribution for an oblique linearly polarized Alfven wave.

    Note: field with a wave propagating along the z axis at a large angle `theta` with respect to the background magnetic field lying in the x-z plane.

    The initial waveis an Alfven mode in which the magnetic field fluctuation points along the y and z axis and has a relative amplitude $A = \delta B_y / B_0$
    """

    Bz_expression = f"{B0 * cosd(theta)}"
    By_expression = f"{A * B0} * cos({k} * z)"
    Bx_expression = f"{B0 * sind(theta)}"
    pz_expression = 0
    # py_expression = f"{A * vA * cosd(theta)} * cos({k} * z)"
    py_expression = f"{A * vA} * cos({k} * z)"
    px_expression = f"{vA * sind(theta)}"

    field = picmi.AnalyticInitialField(
        Bx_expression=Bx_expression,
        By_expression=By_expression,
        Bz_expression=Bz_expression,
    )

    momentum_expressions = [px_expression, py_expression, pz_expression]
    rms_velocity = [v_ti, v_ti, v_ti]
    dist = picmi.AnalyticDistribution(
        density_expression=n0,
        momentum_expressions=momentum_expressions, #: Analytic expressions describing the gamma*velocity for each axis [m/s]
        warpx_momentum_spread_expressions=rms_velocity,  #: gamma*velocity spread for each axis
    )

    return field, dist


class AlfvenModes(HybridSimulation):
    # Applied field parameters
    dim: int = 1
    B0: float = 100 * 1e-9
    """Initial magnetic field strength (T)"""
    n0: float = 100 * 1e6
    """Initial plasma density (m^-3)"""

    A: float = 1  # relative amplitude
    theta: float = 0  # angle with respect to the background magnetic field
    wave_number: int = 1  # wave number

    # Spatial domain
    Lz_norm: float = 128
    nx: int = 16  # by default blocking_factor is 8 so at least 16
    ny: int = 16

    @property
    def _w_ci(self):
        """Modified Ion cyclotron frequency (rad/s) due to large amplitude waves"""
        return constants.q_e * abs(self.B0) * np.sqrt(1+self.A**2) / self.m_ion

    def setup_init_cond(self):
        """setup initial conditions"""
        self.k = 2 * np.pi / (self.Lz / self.wave_number)  # angular wavenumber

        B_ext, dist = init_obl_alfven(
            k=self.k,
            B0=self.B0,
            A=self.A,
            theta=self.theta,
            vA=self.vA,
            n0=self.n0,
            v_ti=self.v_ti,
        )
        self._B_ext = B_ext
        self._dist = dist

        return self


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
    wave_number: float = 1,
    dz_norm: float = 0.5,
    dt_norm: float = 1 / 64,
    time_norm: float = 100,
    substeps: int = 16,
    nppc: int = 64,
    dry_run: bool = False,
):

    base_dir = Path(__file__).parent
    sub_dir = f"dim_{dim}_beta_{beta}_theta_{theta}_eta_{eta}"
    directory = base_dir / sub_dir

    os.makedirs(directory, exist_ok=True)
    os.chdir(directory)

    grid_kwargs = dict()
    if dim == 3:
        grid_kwargs.update(
            warpx_blocking_factor_x=4,
            warpx_blocking_factor_y=4,
        )

    simulation = AlfvenModes(
        **ctx.params,
        plasma_resistivity=eta,
        diag_part=True,
        # Reduce the number of cells in x and y to accelerate the simulation
        grid_kwargs=grid_kwargs,
        nx=8,
        ny=8,
    )

    print(ctx.params)
    # wavenumber
    print("wavenumer (1 / ion inertial length):", simulation.k * simulation.d_i)

    # wave period
    w = simulation.k * simulation.vA
    t_w = 2 * np.pi / w
    t_w_norm = t_w / simulation.t_ci
    t_w_norm

    if not dry_run:
        simulation._sim.step()

    return simulation

if __name__ == "__main__":
    app()
