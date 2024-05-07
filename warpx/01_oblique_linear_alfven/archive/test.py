# %%
import os
import yt
from functools import partial

dim = 1
directory = "dim_1_beta_0.25_theta_60_eta_200"
field_diag_dir = "diags/diag1"
file_format = "??????"

# %%
os.chdir(directory)
ts_field = yt.load(f"{field_diag_dir}{file_format}")
ds_field = ts_field[0]

# because `yt` will always load the data in 3D, we need to specify the direction corresponding to the simulation direction `z`
if dim == 3:
    direction = "z"
elif dim == 2:
    direction = "y"
else:
    direction = "x"

# %%
plot_field = partial(
    yt.ProfilePlot,
    x_field=direction,
    weight_field=("boxlib", "volume"),
    x_log=False,
    y_log=False,
)

plot_field(ds_field, y_fields="By")


