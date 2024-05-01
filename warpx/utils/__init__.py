from yt.data_objects.time_series import SimulationTimeSeries
from yt.data_objects.static_output import Dataset
from yt.frontends.boxlib.data_structures import WarpXDataset
import pandas as pd
import xarray as xr

from xarray import DataArray
from pandas import DataFrame

def rename(data: DataFrame | DataArray | xr.Dataset, dict):
    if isinstance(data, (DataArray, xr.Dataset)):
        return data.rename(dict)
    elif isinstance(data, DataFrame):
        return data.rename(columns=dict)

def rename_coords(data: DataFrame | DataArray, ds: Dataset):
    dim = ds.dimensionality
    if isinstance(ds, WarpXDataset):
        if dim == 1 or dim == 2:
            mapper = {"x": "z", "y": "x", "z": "y", "particle_position_x": "particle_position_z"}
            data = rename(data, mapper)
    return data


def check_ds_type(ds: Dataset):
    if len(ds.particle_types) > 0:
        return "particle"
    else:
        return "field"


def ds2df(ds: Dataset, type=None, coords=None):
    type = type or check_ds_type(ds)
    if type == "particle":
        fields = ds.field_list
    elif type == "field":
        coords = coords or ["x", "y", "z"]
        fields = ds.field_list + coords
    df: pd.DataFrame = ds.all_data().to_dataframe(fields)
    return df.assign(time=ds.current_time).pipe(rename_coords, ds)


def export_ds(ds: Dataset, **kwargs):
    # new_fn = ds.filename(ext, "arrow")
    new_fn = ds.filename + ".arrow"
    df = ds2df(ds, **kwargs)
    df.to_feather(new_fn)
    return new_fn


def export_ts(ts: SimulationTimeSeries, **kwargs):
    return [export_ds(ds, **kwargs) for ds in ts.piter()]


def ds2xr(ds: Dataset, fields=None):
    fields = fields or ds.field_list
    grid = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    da: DataArray = grid.to_xarray(fields)
    return da.assign_coords(time=ds.current_time).pipe(rename_coords, ds)

def ts2xr(ts, fields=None):
    fields = fields or ts[0].field_list
    return xr.concat((ds2xr(ds, fields) for ds in ts), dim="time")
