import yt
from yt.data_objects.time_series import SimulationTimeSeries
from yt.data_objects.static_output import Dataset
from yt.frontends.boxlib.data_structures import WarpXDataset
import pandas as pd
import xarray as xr
from tqdm import tqdm

from xarray import Dataset as xrDataset
from pandas import DataFrame
import logging


ytLogger = logging.getLogger("yt")
ytLogger.setLevel(logging.WARNING)

def load_ts(meta, type="field"):
    diag_format = meta.get("diag_format")
    if diag_format == "openpmd":
        ext = meta.get("diag_openpmd_backend")
        file_format = f"/{diag_format}*.{ext}"
    else:
        file_format = "??????"  # Consider more appropriate fallback or error handling

    paths = {"field": "diags/diag1", "particle": "diags/diag2"}

    return yt.load(f"{paths[type]}{file_format}")


def load_ts_all(meta):
    return load_ts(meta, "field"), load_ts(meta, "particle")


from plum import dispatch


@dispatch
def _rename(data: xrDataset, dict: dict):
    return data.rename(dict)


@dispatch
def _rename(data: DataFrame, dict: dict):
    return data.rename(columns=dict)


def rename_coords(data, ds):
    dim = ds.dimensionality
    if check_ds_type(ds) == "field":
        prefix = ""
    else:
        prefix = "particle_position_"

    if isinstance(ds, WarpXDataset):
        if dim == 1 or dim == 2:
            mapper = {
                f"{prefix}x": f"{prefix}z",
                f"{prefix}y": f"{prefix}x",
                f"{prefix}z": f"{prefix}y",
            }
            data = _rename(data, mapper)
    return data


def check_ds_type(ds: Dataset):
    if len(ds.particle_types) > 0:
        return "particle"
    else:
        return "field"


def ds2df(ds: Dataset, fields: list = None, coords=["x", "y", "z"]):
    fields = fields or ds.field_list
    if check_ds_type(ds) == "field":
        fields = fields + coords
    df: pd.DataFrame = ds.all_data().to_dataframe(fields)
    return df.assign(time=ds.current_time).pipe(rename_coords, ds)


def export_ds(ds: Dataset, **kwargs):
    # new_fn = ds.filename(ext, "arrow")
    new_fn = ds.filename + ".arrow"
    df = ds2df(ds, **kwargs)
    df.to_feather(new_fn)
    return new_fn


def export_ts(ts: SimulationTimeSeries, **kwargs):
    return [export_ds(ds, **kwargs) for ds in tqdm(ts.piter())]


def ds2xr(ds: Dataset, fields: list = None) -> xr.Dataset:
    fields = fields or ds.field_list
    grid = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    da: xr.Dataset = grid.to_xarray(fields)
    return da.assign_coords(time=ds.current_time).pipe(rename_coords, ds)


def ts2xr(ts, fields=None):
    fields = fields or ts[0].field_list
    return xr.concat((ds2xr(ds, fields) for ds in ts), dim="time")
