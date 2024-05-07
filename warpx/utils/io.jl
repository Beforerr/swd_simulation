using DataFrames,
    DataFramesMeta
using Arrow
using PhysicalConstants.CODATA2018: ElementaryCharge
using Unitful

function unit_df!(df)
    @chain df begin
        @transform!(
            :rho_c = :rho .* 1u"C/m^3",
            :velocity_th_parp = :velocity_th_parp .* 1u"m/s",
            :velocity_th_perp = :velocity_th_perp .* 1u"m/s",
        )
        @transform!(
            :rho_n = :rho_c ./ ElementaryCharge,

            # T_plasma
        )
    end 
end


function normalize_df!(df, meta)
    n0 = meta["n0"] * u"1/m^3"
    T0 = meta["T_plasma"] * u"eV"
    mass = meta["m_ion"] * u"kg"
    @transform!(df,
        :rho_n_norm = :rho_n ./ n0,
        :time_norm = :time ./ meta["t_ci"],
        :z_norm = :z ./ meta["d_i"],
        :T_parp_norm = mass * :velocity_th_parp.^2 / T0,
        :T_perp_norm = mass * :velocity_th_perp.^2 / T0,
    )
end

function process_df!(df)
    transform!(df,
        names(df, r"B") => ByRow(norm ∘ vcat) => :Bmag,
        names(df, r"E") => ByRow(norm ∘ vcat) => :Emag,
        names(df, r"j") => ByRow(norm ∘ vcat) => :jmag,
    )
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
    normalize_df!(df, meta) |> process_df!
    println(names(df))
    sort!(df, [:time, :z, :y, :x])
end