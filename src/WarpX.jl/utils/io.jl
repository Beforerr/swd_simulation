using DataFrames,
    DataFramesMeta,
    CategoricalArrays
using FileIO
using Arrow
import JSON
using PhysicalConstants.CODATA2018: ElementaryCharge, μ_0
# using Statistics
using LinearAlgebra
using Unitful

include("warpx.jl")

"""
load simulation metadata (json)
"""
function load_meta(filename="sim_parameters.json")
    # data = JSON.parsefile(filename)
    data = load(filename)
    return data
end

"""
load simulation output field data
"""
function load_output_field(meta)
    field_diag_dir = (meta["diag_format"] == "openpmd" ? "diags/diag1" : "diags")
    files = filter(contains(r".*\.arrow"), readdir(field_diag_dir, join=true))
    dfs = files .|> Arrow.Table .|> DataFrame
    reduce(vcat, dfs)
end

velocity_comps = [:velocity_th_x, :velocity_th_y, :velocity_th_z, :velocity_th_parp, :velocity_th_perp]
T_comps = [:T_x, :T_y, :T_z, :T_parp, :T_perp]
T_norm_comps = [:T_x_norm, :T_y_norm, :T_z_norm, :T_parp_norm, :T_perp_norm]

function unit_df!(df)
    @chain df begin
        transform!(
            velocity_comps .=> v -> v .* 1u"m/s", renamecols=false
        )
        @transform!(
            :rho_c = :rho .* 1u"C/m^3",
        )
        @transform!(
            :rho_n = :rho_c ./ ElementaryCharge,
        )
    end
end


function transform_pair(func, x, y)
    return x => func => y
end

function transform_map(func, x, y)
    return map(transform_pair$func, x, y)
end

function normalize_df!(df, meta)
    n0 = meta["n0"] * u"1/m^3"
    T0 = meta["T_plasma"] * u"eV"
    mass = meta["m_ion"] * u"kg"

    function normalize_temp(x)
        convert.(Float64, mass * x .^ 2 / T0)
    end

    @chain df begin
        @transform!(
            :rho_n_norm = :rho_n ./ n0,
            :time_norm = :time ./ meta["t_ci"],
            :z_norm = :z ./ meta["d_i"],
        )
        transform!(
            transform_map(normalize_temp, velocity_comps, T_norm_comps)...
        )
    end
end

function process_df!(df, meta)
    mass = meta["m_ion"] * u"kg"

    function vth2temp(x)
        mass * x .^ 2
    end

    @chain df begin
        transform!(
            B_comps => ByRow(norm ∘ vcat) => :Bmag,
            E_comps => ByRow(norm ∘ vcat) => :Emag,
            j_comps => ByRow(norm ∘ vcat) => :jMag,
            # j_e_comps => ByRow(norm ∘ vcat) => :jeMag,
            transform_map(vth2temp, velocity_comps, T_comps)...
        )
        @transform!(
            :vA_x = :Bx * u"T" ./ sqrt.(μ_0 * mass * :rho_n),
            :vA_y = :By * u"T" ./ sqrt.(μ_0 * mass * :rho_n),
            :vA_z = :Bz * u"T" ./ sqrt.(μ_0 * mass * :rho_n),
            :Λ_temp = :T_perp ./ :T_parp,
            :Λ = convert.(Float64, μ_0 * :rho_n .* (:T_parp .- :T_perp) ./ (:Bmag * u"T") .^ 2),
        )
        @transform!(
            :Λ_temp_log = log.(:Λ_temp),
            # pressure anisotropy modified Alfven speed
            :vA_p_x = :vA_x .* sqrt.(1 .- :Λ),
            :vA_p_y = :vA_y .* sqrt.(1 .- :Λ),
            :vA_p_z = :vA_z .* sqrt.(1 .- :Λ),
        )
    end
end

function load_pressure_df(filename="pressure.arrow")
    # load(filename) |> DataFrame
    Arrow.Table(filename) |> DataFrame
end

function load_field(meta::AbstractDict)
    df = load_output_field(meta)
    try
        df_p = load_pressure_df()
        df = innerjoin(df, df_p, on=["time", "x", "y", "z"], makeunique=true)
    catch
    end
    unit_df!(df)
    process_df!(df, meta)
    normalize_df!(df, meta)
    sort!(df, [:time, :z, :y, :x])
end

load_field() = load_field(load_meta())
load_field(path::AbstractString) = cd(load_field, path)


using DrWatson

function collect_results(folder=datadir();
    valid_filetypes=[".json"],
    subfolders=true,
    white_list=["dim", "beta", "theta", "plasma_resistivity", "wave_length", "path"]
)

    results = DrWatson.collect_results(folder, subfolders=subfolders, valid_filetypes=valid_filetypes)
    return @chain results begin
        select!(white_list)
        @transform!(
            :df = load_field.(dirname.(:path))
            # :dir = dirname.(:path),
            # :meta = load_meta.(:path),
        )
    end
end
