import numpy as np
import unyt as u
from .io import check_ds_type
from yt.data_objects.static_output import Dataset
from unyt.physical_constants import mu_0

def add_V_Alfven_field(ds: Dataset, direction: str, ion_mass: float = 1.0):
    #: NOTE: rho is species charge density

    def _V_Alfven(field, data):
        rho_m = data["rho"] * ion_mass * u.m ** (-3) / u.qp_mks.value * u.kg
        return data[f"B{direction}"] / np.sqrt(mu_0 * rho_m)

    name = ("boxlib", f"V_Alfven_{direction}")
    ds.add_field(name, function=_V_Alfven, units="km/s", sampling_type="local")


def add_V_field(ds: Dataset, direction: str, ion_mass: float = 1.0):
    def _V(field, data):
        return data[f"particle_momentum_{direction}"] / ion_mass / u.kg

    name = ("ions", f"V_{direction}")
    ds.add_field(name, function=_V, units="km/s", sampling_type="particle")


def add_field(ds: Dataset, meta: dict):
    type = check_ds_type(ds)
    if type == "particle":
        add_V_field(ds, "x", meta["m_ion"])
        add_V_field(ds, "y", meta["m_ion"])
        add_V_field(ds, "z", meta["m_ion"])
    elif type == "field":
        add_V_Alfven_field(ds, "x", meta["m_ion"])
        add_V_Alfven_field(ds, "y", meta["m_ion"])
        add_V_Alfven_field(ds, "z", meta["m_ion"])
