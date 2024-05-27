using DrWatson
@quickactivate "IDsHybridSimulation"

using AlgebraOfGraphics,
    CairoMakie
using DataFrames,
    DataFramesMeta,
    CategoricalArrays
import JSON
using Statistics
using LinearAlgebra
using beforerr
using LaTeXStrings
using PartialFunctions
using Comonicon

include("utils/io.jl")
include("utils/analysis.jl")
include("utils/plot.jl")

# set_aog_theme!()

function setup(directory; base_dir=datadir())
    # change to simulation directory and load metadata
    try
        cd(joinpath(base_dir, directory))
        println("Changed to $directory")
    catch
    end
    # load simulation metadata (json)
    JSON.parsefile("sim_parameters.json")
end

"""
# Options

- `--opt1 <arg>`: an option
- `-o, --opt2 <arg>`: another option

# Flags
"""
@cast function run(;
    directory::String=".",
    noPlotFields::Bool=false, noPlotFieldsTime::Bool=false, noPlotOverviewTs::Bool=false
    # dim::Int=1, beta::Float64=0.25, theta::Float64=60.0, eta::Float64=10.0,
)
    meta = setup(directory)
    df = load_field(meta)

    noPlotFields || plot_fields(df)
    noPlotFieldsTime || plot_fields_time(df)
    noPlotOverviewTs || plot_overview_ts(df)
end


if abspath(PROGRAM_FILE) == @__FILE__
    @main
end


# # calculate the mean of the data by averaging over "y" and "z"

# function plot_fields(df, fields; ids=ids, fig_options = (size = (800, 800),))

#     temp_df = get_avg_fields(df, fields, ids=ids)
#     plt = data(temp_df) * mapping(Pair.(ids, labs)..., :value, row = :variable) * visual(Heatmap)

#     draw(plt; figure = fig_options)
# end

# plot_fields(df, B_fields)
# easy_save("B_field")
# plot_fields(df, E_field)
# easy_save("E_field")
# plot_fields(df, j_field)
# easy_save("j_field")
# plot_fields(df, rho_field)
# easy_save("rho_field");


# pressure_cols = names(df, r"pressure")
# temp_cols = replace.(pressure_cols, r"pressure" => "T")