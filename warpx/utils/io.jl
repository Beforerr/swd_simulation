using DataFrames,
    DataFramesMeta
using Arrow
using PhysicalConstants.CODATA2018: ElementaryCharge, μ_0
using Unitful

B_components = ["Bx", "By", "Bz"]
E_components = ["Ex", "Ey", "Ez"]
j_components = ["jx", "jy", "jz"]
velocity_components = [:velocity_th_x, :velocity_th_y, :velocity_th_z, :velocity_th_parp, :velocity_th_perp]
T_components = [:T_x, :T_y, :T_z, :T_parp, :T_perp]
T_norm_components = [:T_x_norm, :T_y_norm, :T_z_norm, :T_parp_norm, :T_perp_norm]

function unit_df!(df)
    @chain df begin
        transform!(
            velocity_components .=> v -> v .* 1u"m/s", renamecols=false
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
            transform_map(normalize_temp, velocity_components, T_norm_components)...
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
            B_components => ByRow(norm ∘ vcat) => :Bmag,
            E_components => ByRow(norm ∘ vcat) => :Emag,
            j_components => ByRow(norm ∘ vcat) => :jmag,
            transform_map(vth2temp, velocity_components, T_components)...
        )
        @transform!(
            :vA_x = :Bx * u"T" ./ sqrt.(μ_0 * mass * :rho_n),
            :vA_y = :By * u"T" ./ sqrt.(μ_0 * mass * :rho_n),
            :vA_z = :Bz * u"T" ./ sqrt.(μ_0 * mass * :rho_n),
            :Λ_temp = :T_perp ./ :T_parp,
            :Λ = convert.(Float64, μ_0 * :rho_n .* (:T_parp .- :T_perp ) ./ (:Bmag * u"T").^2 ),
        )
        @transform!( # pressure anisotropy modified Alfven speed
            :vA_p_x = :vA_x .* sqrt.(1 .- :Λ),
            :vA_p_y = :vA_y .* sqrt.(1 .- :Λ),
            :vA_p_z = :vA_z .* sqrt.(1 .- :Λ),
        )
    end
end

function load_output_field(meta)
    field_diag_dir = (meta["diag_format"] == "openpmd" ? "diags/diag1" : "diags")
    files = filter(contains(r".*\.arrow"), readdir(field_diag_dir, join=true))
    dfs = files .|> Arrow.Table .|> DataFrame
    reduce(vcat, dfs)
end

function load_pressure_df(fp="pressure.arrow")
    Arrow.Table(fp) |> DataFrame
end

function load_field(meta;)
    df = load_output_field(meta)
    try
        df_p = load_pressure_df()
        df = innerjoin(df, df_p, on=["time", "x", "y", "z"], makeunique=true)
    catch
    end
    unit_df!(df)
    process_df!(df, meta)
    normalize_df!(df, meta)
    println(names(df))
    sort!(df, [:time, :z, :y, :x])
end