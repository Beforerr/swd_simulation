# %% [markdown]
# # Analysis

# %%
using AlgebraOfGraphics,
    CairoMakie
using DataFrames,
    DataFramesMeta,
    CategoricalArrays
using Arrow
using Statistics
import JSON

# %%
dim = 3
beta = 0.25
theta = 60
dir = "$(@__DIR__)/01_oblique_linear_alfven/dim_$(dim)_beta_$(beta)_$(theta)_60"
# load simulation metadata (json)
meta = JSON.parsefile("$(dir)/sim_parameters.json")

println("Simulation parameters:")
for (key, value) in meta
    println("  $key: $value")
end

# %% [markdown]
# # Particle data

# %%
file = "data.arrow"
path = joinpath(dir, file)
df = path |> Arrow.Table |> DataFrame

# %%
B_fields = names(df, r"B")
E_field = names(df, r"E")
j_field = names(df, r"j")

variables = [B_fields; E_field; j_field]

# %%
# calculate the mean of the data by averaging over "y" and "z"
ids = [:z, :time]

function plot_fields(df, fields; ids=ids)
    temp_df = @chain df begin
        groupby(ids)
        combine(fields .=> mean, renamecols=false)
        stack(fields, ids)
    end

    plt = data(temp_df) * mapping(ids..., :value, row = :variable) * visual(Heatmap)
    draw(plt)
end

# %%
plot_fields(df, B_fields)

# %%
plot_fields(df, E_field)

# %%
plot_fields(df, j_field)

# %% [markdown]
# ## Fluid fields

# %%
file = "particle.arrow"
path = joinpath(dir, file)
df = path |> Arrow.Table |> DataFrame

# %%
df.z_norm = df.particle_position_z / meta["d_i"]
df.py_norm = df.particle_momentum_y / 1e-25

# %%
plt = data(df) * mapping(:particle_position_z, :time,:particle_momentum_y) * visual(Heatmap)
draw(plt)

# %% [markdown]
# ### Non-binned particle data

# %%
using CategoricalArrays

# %%
df.time_norm = CategoricalArray(df.time ./ meta["t_ci"])

# %%
z_norm_edge = 0:1:240
py_norm_edge = -1e3:10:1e3

# %%
datalimits_f = x -> quantile(x, [0.05, 0.95])

# %%
fig_options = (size = (1200, 1000),)

plt = data(df) * mapping(:z_norm, :py_norm, layout=:time_norm) * histogram(datalimits=datalimits_f)
p = draw(plt; figure = fig_options)

# %% [markdown]
# ## Parameters

# %%
using Pkg
Pkg.add("Symbolics")

# %%
using Symbolics

# %%
@variables t x y μ_0 B ρ c n q m ϵ_0

Alfven_speed = B / sqrt(μ_0 * ρ)

# plasma frequency
ω_p = sqrt(n * q / (m * ϵ_0))
# inertial_length
d_i = c / ω_p

# gryofrequency
ω_c = q * B / m


simplify(ω_c / ω_p)

# %%
ω_c / ω_p

