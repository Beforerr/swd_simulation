import numpy as np
from pywarpx import picmi
from space_analysis.utils.math import cosd, sind
from space_analysis.simulation.warpx import HybridSimulation

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
        momentum_expressions=momentum_expressions,  #: Analytic expressions describing the gamma*velocity for each axis [m/s]
        warpx_momentum_spread_expressions=rms_velocity,  #: gamma*velocity spread for each axis
    )

    return field, dist


class AlfvenModes(HybridSimulation):
    # Applied field parameters
    dim: int = 1
    n0: float = 100 * 1e6
    """Initial plasma density (m^-3)"""

    A: float = 1  # relative amplitude
    theta: float = 0  # angle with respect to the background magnetic field
    wave_length: float  # wave length normalized to the ion inertial length

    # Spatial domain
    Lz_norm: float = 64
    nx: int = 16  # by default blocking_factor is 8 so at least 16
    ny: int = 16

    @property
    def _w_ci(self):
        """Modified Ion cyclotron frequency (rad/s) due to large amplitude waves"""
        return constants.q_e * abs(self.B0) * np.sqrt(1 + self.A**2) / self.m_ion

    def setup_init_cond(self):
        """setup initial conditions"""
        self.k = 2 * np.pi / (self.wave_length * self.d_i)  # angular wavenumber

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

