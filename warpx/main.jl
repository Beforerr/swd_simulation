using AlgebraOfGraphics,
    CairoMakie
using DataFrames,
    DataFramesMeta,
    CategoricalArrays
using Arrow
using Statistics
using LinearAlgebra
using LaTeXStrings
using Latexify
using beforerr
using PartialFunctions
include("utils/io.jl")
include("utils/analysis.jl")

dim = 3
dim = 1
beta = 0.25
theta = 60.0
eta = 100.0
eta = 10.0

meta = setup(dim, beta, theta, eta)
df = load_field(meta);

B_fields = names(df, r"^B(x|y|z|mag)$")
E_field = names(df, r"E")
j_field = names(df, r"j")
rho_field = names(df, r"rho")

variables = [B_fields; E_field; j_field]


# calculate the mean of the data by averaging over "y" and "z"
z_norm_lab = L"z~(d_i)"
time_norm_lab = L"T~(T_{ci})"
ids = [:z_norm, :time_norm]
labs = [z_norm_lab, time_norm_lab]
function plot_fields(df, fields; ids=ids, fig_options = (size = (800, 800),))

    temp_df = get_avg_fields(df, fields, ids=ids)
    plt = data(temp_df) * mapping(Pair.(ids, labs)..., :value, row = :variable) * visual(Heatmap)

    draw(plt; figure = fig_options)
end

plot_fields(df, B_fields)
easy_save("B_field")
plot_fields(df, E_field)
easy_save("E_field")
plot_fields(df, j_field)
easy_save("j_field")
plot_fields(df, rho_field)
easy_save("rho_field");


pressure_cols = names(df, r"pressure")
temp_cols = replace.(pressure_cols, r"pressure" => "T")


add_labels!