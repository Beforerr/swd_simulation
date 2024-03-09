#!/usr/bin/env python
# coding: utf-8

# In[12]:


import os

import numpy as np
from mpi4py import MPI as mpi
from pywarpx import callbacks, fields, picmi
from pywarpx import libwarpx
from pywarpx.picmi import Simulation

import dill

constants = picmi.constants

comm = mpi.COMM_WORLD


# In[13]:


B0 = 1.  # Initial magnetic field strength (T)
pert = 1  # Perturbation amplitude


# In[14]:


def dump(x, file="sim_parameters.dpkl"):
    """dump all the current attributes to a dill pickle file"""
    if comm.rank == 0:
        with open(file, "wb") as f:
            dill.dump(x, f)


# In[15]:


class CustomSimulation(Simulation):

    dim: int = None

    # Plasma parameters
    n0 = None
    """plasma density (m^-3)"""

    # Plasma species parameters
    m_ion_norm: float = None
    """Ion mass (electron masses)"""
    m_ion: float = None
    """Ion mass (kg)"""
    v_ti: float = None
    """Ion thermal velocity (m/s)"""
    

    # Spatial domain
    nz: int = None  #: number of cells in z direction
    nx: int = None  #: number of cells in x (and y) direction for >1 dimensions

    # Numerical parameters
    nppc: int = None
    """Seed number of particles per cell"""
    diag_steps: int = None
    """The simulation step interval at which to output diagnostic"""

    def post_init(self):
        """This function is called after the object is initialized"""
        pass

    def setup_particle(self):
        """setup the particle"""
        self.ions = picmi.Species(
            name="ions",
            charge_state=1,
            mass=self.m_ion,
            initial_distribution=picmi.UniformDistribution(
                density=self.n0,
                rms_velocity=[self.v_ti] * 3,
            ),
        )
        self.add_species(
            self.ions,
            layout=picmi.PseudoRandomLayout(
                grid=self.grid, n_macroparticles_per_cell=self.nppc
            ),
        )
        return self

    def setup_grid(self):
        """Setup geometry and boundary conditions"""
        if self.dim == 1:
            grid_object = picmi.Cartesian1DGrid
        elif self.dim == 2:
            grid_object = picmi.Cartesian2DGrid
        else:
            grid_object = picmi.Cartesian3DGrid

        number_of_cells = [self.nx, self.nx, self.nz][-self.dim :]
        boundary_conditions = ["periodic"] * self.dim

        self.grid = grid_object(
            number_of_cells=number_of_cells,
            warpx_max_grid_size=self.nz,
            lower_bound=[-self.Lx / 2.0, -self.Lx / 2.0, 0][-self.dim :],
            upper_bound=[self.Lx / 2.0, self.Lx / 2.0, self.Lz][-self.dim :],
            lower_boundary_conditions=boundary_conditions,
            upper_boundary_conditions=boundary_conditions,
        )
        return self    

    def setup_field(self):
        pass

    def setup_diag(self):
        self.fmt = "txt"
        self.diag_steps = int(1 / 4 / self.dt_norm)
        
        
        if self.B_dir == "z":
            self.output_file_name = "par_field_data"
        else:
            self.output_file_name = "perp_field_data"

        if self.B_dir == "z" or self.dim == 1:
            line_diag = picmi.ReducedDiagnostic(
                diag_type="FieldProbe",
                probe_geometry="Line",
                z_probe=0,
                z1_probe=self.Lz,
                resolution=self.nz - 1,
                name=self.output_file_name,
                period=self.diag_steps,
                path="diags/",
            )
            self.add_diagnostic(line_diag)
        else:
            # install a custom "reduced diagnostic" to save the average field
            callbacks.installafterEsolve(self._record_average_fields)
            try:
                os.mkdir("diags")
            except OSError:
                # diags directory already exists
                pass
            with open(f"diags/{self.output_file_name}.{self.fmt}", "w") as f:
                f.write(
                    "[0]step() [1]time(s) [2]z_coord(m) "
                    "[3]Ez_lev0-(V/m) [4]Bx_lev0-(T) [5]By_lev0-(T)\n"
                )

        field_diag = picmi.FieldDiagnostic(
            grid=self.grid,
            period=self.diag_steps,
            # warpx_format="openpmd",
            # warpx_openpmd_backend = 'h5'
        )
        self.add_diagnostic(field_diag)

    def setup_run(self):
        """Setup simulation components."""
        self.setup_grid().setup_field()
        self.setup_particle()
        self.setup_diag()
        self.initialize_inputs()

    def _record_average_fields(self):
        """A custom reduced diagnostic to store the average E&M fields in a
        similar format as the reduced diagnostic so that the same analysis
        script can be used regardless of the simulation dimension.
        """
        step = self.extension.warpx.getistep(lev=0) - 1

        # Skip steps not aligned with diagnostic steps
        if step % self.diag_steps != 0:
            return

        if libwarpx.amr.ParallelDescriptor.MyProc() != 0:
            return

        Bx_warpx = fields.BxWrapper()[...]
        By_warpx = fields.ByWrapper()[...]
        Ez_warpx = fields.EzWrapper()[...]

        t = step * self.time_step_size
        z_vals = np.linspace(0, self.Lz, self.nz, endpoint=False)

        if self.dim == 2:
            Ez = np.mean(Ez_warpx[:-1], axis=0)
            Bx = np.mean(Bx_warpx[:-1], axis=0)
            By = np.mean(By_warpx[:-1], axis=0)
        else:
            Ez = np.mean(Ez_warpx[:-1, :-1], axis=(0, 1))
            Bx = np.mean(Bx_warpx[:-1], axis=(0, 1))
            By = np.mean(By_warpx[:-1], axis=(0, 1))

        # Prepare output lines
        output_lines = [
            f"{step:05d} {t:.10e} {z:.10e} {Ez_val:+.10e} {Bx_val:+.10e} {By_val:+.10e}\n"
            for z, Ez_val, Bx_val, By_val in zip(z_vals, Ez, Bx, By)
        ]

        with open(f"diags/{self.output_file_name}.{self.fmt}", "a") as f:
            f.writelines(output_lines)


# In[16]:


class HybridSimulation(CustomSimulation):

    beta: float = None
    """Plasma beta"""  # used to calculate temperature
    
    vA: float = None  # Alfven speed
    vA_over_c: float = None  # Alfven speed over speed of light

    plasma_resistivity: float = None  # Plasma resistivity

    T_plasma: float = None
    Te: float = None
    """Electron temperature in (eV)"""

    t_ci: float = None
    """Ion cyclotron period (s)"""
    l_i: float = None
    """Ion inertial length (m)"""


    # Numerical parameters
    time_norm = 300.0
    """Simulation temporal length (ion cyclotron periods)"""
    dt_norm: float = 1 / 128
    """Time step (ion cyclotron periods)"""

    dz_norm: float = 1 / 16
    """Cell size (ion skin depths)"""
    nz: int = None
    """number of cells in z direction"""

    substeps = None  # Number of substeps used to update B

    def post_init(self):
        """This function is called after the object is initialized"""
        super().post_init()      
        self.time_step_size = self.dt_norm * self.t_ci
        
        self.dz = self.dz_norm * self.l_i
        self.Lz = self.nz * self.dz
        self.Lx = self.nx * self.dz

    def setup_field(self):
        """Setup external field and field solver"""

        self.solver = picmi.HybridPICSolver(
            grid=self.grid,
            Te=self.Te,
            n0=self.n0,
            plasma_resistivity=self.plasma_resistivity,
            substeps=self.substeps,
        )
        
        k1 = 2 * np.pi / self.Lz
        
        B_ext = picmi.AnalyticInitialField(
            Bx_expression="pert*(sin(k1*z + phi1))",
            By_expression="pert*(cos(k1*z + phi1))",
            Bz_expression=B0,
            pert=pert,
            k1=k1,
            phi1=0,
            # By_expression="pert(sin(k1*x + phi1) + sin(k2*x + phi2) + sin(k3*x + phi3))",
            # k2 = 2*np.pi/0.02,
            # phi2 = 0,
        )

        self.add_applied_field(B_ext)
        return self

    def setup_run(self):
        self.current_deposition_algo = "direct"
        self.particle_shape = 1
        super().setup_run()


# In[17]:


def log_info(sim: HybridSimulation):
    """print out plasma parameters and numerical parameters."""
    if comm.rank == 0:
        print(
            f"Initializing simulation with input parameters:\n"
            f"\tTe = {sim.Te:.3f} eV\n"
            f"\tn = {sim.n0:.1e} m^-3\n"
            f"\tB0 = {sim.B0:.2f} T\n"
            f"\tM/m = {sim.m_ion_norm:.0f}\n"
        )
        print(
            f"Plasma parameters:\n"
            f"\tl_i = {sim.l_i:.1e} m\n"
            f"\tt_ci = {sim.t_ci:.1e} s\n"
            f"\tv_ti = {sim.v_ti:.1e} m/s\n"
            f"\tvA = {sim.vA:.1e} m/s\n"
        )
        print(
            f"Numerical parameters:\n"
            f"\tdz = {sim.dz:.1e} m\n"
            f"\tdt = {sim.time_step_size:.1e} s\n"
            f"\tdiag steps = {sim.diag_steps:d}\n"
            f"\ttotal steps = {sim.max_steps:d}\n"
        )


# In[18]:


class EMModes(HybridSimulation):
    """The following runs a simulation of an uniform plasma at a set
    temperature (Te = Ti) with an external magnetic field applied in either the
    z-direction (parallel to domain) or x-direction (perpendicular to domain).
    The analysis script (in this same directory) analyzes the output field data
    for EM modes. This input is based on the EM modes tests as described by
    Munoz et al. (2018) and tests done by Scott Nicks at TAE Technologies.
    """

    # Applied field parameters
    B0 = B0
    betas = [0.01, 0.1]  # Plasma beta, used to calculate temperature

    # Spatial domain
    nx: int = 8  # number of cells in x (and y) direction for >1 dimensions

    # Plasma resistivity - used to dampen the mode excitation
    eta = [[1e-7, 1e-7], [1e-7, 1e-5], [1e-7, 1e-4]]

    l_i: float = None
    """ion skin depth"""

    B_dir: str = None

    _m_ion_norms = [100.0, 400.0]
    _vA_over_cs = [1e-4, 1e-3]
    _nzs = [1024, 1920]  # number of cells in z direction
    # _nppcs = [1024, 256, 64]
    _nppcs = [256, 256, 64]

    def __init__(self, test, dim, B_dir, verbose):
        """Get input parameters for the specific case desired."""
        super().__init__(warpx_serialize_initial_conditions=True)
        self.test = test
        self.dim = dim
        self.B_dir = B_dir
        self.verbose = verbose or self.test

        # get simulation parameters from the defaults given the direction of
        # the initial B-field and the dimensionality
        self.get_simulation_parameters()

        # calculate various plasma parameters based on the simulation input
        self.get_plasma_quantities()

        self.post_init()

        if not self.test:
            self.max_steps = int(self.time_norm / self.dt_norm)
        else:
            # if this is a test case run for only a small number of steps
            self.max_steps = 250

        self.setup_run()
        dump(self)
        log_info(self)

    def get_plasma_quantities(self):
        """Calculate various plasma parameters based on the simulation input."""
        # Ion mass (kg)
        self.m_ion = self.m_ion_norm * constants.m_e

        # Cyclotron angular frequency (rad/s) and period (s)
        self.w_ci = constants.q_e * abs(self.B0) / self.m_ion
        self.t_ci = 2.0 * np.pi / self.w_ci

        # Alfven speed (m/s)
        # vA = B / sqrt(mu0 * n * (M + m)) = c * omega_ci / w_pi
        self.vA = self.vA_over_c * constants.c
        self.n0 = (self.B0 / self.vA) ** 2 / (
            constants.mu0 * (self.m_ion + constants.m_e)
        )

        # Ion plasma frequency (rad/s)
        self.w_pi = np.sqrt(constants.q_e**2 * self.n0 / (self.m_ion * constants.ep0))

        # Skin depth (m): inertial length
        self.l_i = constants.c / self.w_pi

        # Ion thermal velocity (m/s) from beta = 2 * (v_ti / vA)**2
        self.v_ti = np.sqrt(self.beta / 2.0) * self.vA

        # Temperature (eV) from thermal speed: v_ti = sqrt(kT / M)
        self.T_plasma = self.v_ti**2 * self.m_ion / constants.q_e  # eV
        self.Te = self.T_plasma

        # Larmor radius (m)
        self.rho_i = self.v_ti / self.w_ci

    def get_simulation_parameters(self):
        """Pick appropriate parameters from the defaults given the direction
        of the B-field and the simulation dimensionality."""
        if self.B_dir == "z":
            idx = 0
            self.Bx = 0.0
            self.By = 0.0
            self.Bz = self.B0
        elif self.B_dir == "y":
            idx = 1
            self.Bx = 0.0
            self.By = self.B0
            self.Bz = 0.0
        else:
            idx = 1
            self.Bx = self.B0
            self.By = 0.0
            self.Bz = 0.0

        self.beta = self.betas[idx]
        self.m_ion_norm = self._m_ion_norms[idx]
        self.vA_over_c = self._vA_over_cs[idx]
        self.nz = self._nzs[idx]

        self.nppc = self._nppcs[self.dim - 1]
        self.plasma_resistivity = self.eta[self.dim - 1][idx]


# In[19]:


import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument(
    "-t",
    "--test",
    help="toggle whether this script is run as a short CI test",
    action="store_true",
)
parser.add_argument(
    "-d", "--dim", help="Simulation dimension", required=False, type=int, default=1
)
parser.add_argument(
    "--bdir",
    help="Direction of the B-field",
    required=False,
    choices=["x", "y", "z"],
    default="z",
)
parser.add_argument(
    "-v",
    "--verbose",
    help="Verbose output",
    action="store_true",
)
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1] + left

# sanity check
dim = args.dim
assert dim > 0 and dim < 4, f"{dim}-dimensions not a valid input"


# In[20]:


simulation = EMModes(test=args.test, dim=args.dim, B_dir=args.bdir, verbose=args.verbose)


# In[23]:


6/simulation.Lz*0.001


# In[10]:


simulation.write_input_file()
simulation.initialize_warpx()


# In[10]:


simulation.step()


# ## Test

# In[9]:


from plasmapy.formulary import plasma_frequency, inertial_length
import astropy.units as u


# In[1]:


def test():
    plasma_frequency(simulation.n0 * u.m**-3, 'p+')


# In[ ]:




