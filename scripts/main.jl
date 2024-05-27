using DrWatson
@quickactivate "IDsHybridSimulation"

using Statistics
using LinearAlgebra
using beforerr
using LaTeXStrings
using PartialFunctions
using Comonicon

include("utils/io.jl")
include("utils/analysis.jl")
include("utils/plot.jl")

"""
# Options

- `--opt1 <arg>`: an option
- `-o, --opt2 <arg>`: another option

# Flags
"""
@cast function run(;
    directory::String=".",
    noPlotFields::Bool=false, noPlotFieldsTime::Bool=false, noPlotOverviewTs::Bool=false,
    data_dir=datadir()
    # dim::Int=1, beta::Float64=0.25, theta::Float64=60.0, eta::Float64=10.0,
)
    prefix, parameters, suffix = parse_savename(directory)
    base_dir = joinpath(data_dir, directory)
    df = cd(base_dir) do
        load_field()
    end

    noPlotFields || plot_fields(df)
    noPlotFieldsTime || plot_fields_time(df)
    noPlotOverviewTs || plot_overview_ts(df)
end


if abspath(PROGRAM_FILE) == @__FILE__
    @main
end