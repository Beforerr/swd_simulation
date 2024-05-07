using DataFrames,
    DataFramesMeta
using Arrow
import JSON
using PhysicalConstants.CODATA2018: ElementaryCharge
using Unitful

function setup(dim, beta, theta, eta)
    # change to simulation directory and load metadata
    dir = "01_oblique_linear_alfven/dim_$(dim)_beta_$(beta)_theta_$(theta)_eta_$(eta)"
    try
        cd(dir)
    catch
    end
    # load simulation metadata (json)
    JSON.parsefile("sim_parameters.json")
end


function unit_df!(df)
    @chain df begin
        @transform!(
            :rho_c = :rho .* 1u"C/m^3",
        )
        @transform!(
            :rho_n = :rho_c ./ ElementaryCharge,
        )
        @transform!(
            :pressure_perp_u = :pressure_perp .* 1u"kg*(m/s)^2" ./ :rho_n,
            :pressure_parp_u = :pressure_parp .* 1u"kg*(m/s)^2" ./ :rho_n,
        )
    end 
end


function normalize_df!(df, meta)
    n0 = meta["n0"] * u"1/m^3"
    @transform!(df,
        :rho_n_norm = :rho_n ./ n0,
        :time_norm = :time ./ meta["t_ci"],
        :z_norm = :z ./ meta["d_i"],
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
    println(names(df))
    unit_df!(df)
    normalize_df!(df, meta) |> process_df!
    sort!(df, [:time, :z, :y, :x])
end