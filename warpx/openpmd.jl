# %% [markdown]
# # Analysis

# %%
using HDF5

# %%
dim = 3
dir = "$(@__DIR__)/01_oblique_linear_alfven/dim_$(dim)_beta_0.25_theta_60"

filename = joinpath(dir, "diags/diag2/openpmd_000000.h5")
fid = h5open(filename)
data = fid["data"]