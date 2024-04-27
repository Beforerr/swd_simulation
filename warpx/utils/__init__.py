from yt.data_objects.time_series import SimulationTimeSeries
from yt.data_objects.static_output import Dataset
import pandas as pd

def ds2df(ds: Dataset, type = "particle", coords = None) -> pd.DataFrame:
    if type == "particle":
        fields = ds.field_list
    elif type == "field":
        coords = coords or ["x", "y", "z"]
        fields = ds.field_list + coords
    return ds.all_data().to_dataframe(fields).assign(time = ds.current_time)

def export_ds(ds: Dataset, **kwargs):
    # new_fn = ds.filename(ext, "arrow")
    new_fn = ds.filename + ".arrow"
    df = ds2df(ds, **kwargs)
    df.to_feather(new_fn)
    return new_fn

def export_ts(ts: SimulationTimeSeries,  **kwargs):
    return [export_ds(ds, **kwargs) for ds in ts.piter()]