using DataFrames,
    DataFramesMeta
using Arrow
import JSON

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

function normalize_df!(df)
    @transform!(df,
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
    dfs = vcat(files .|> Arrow.Table .|> DataFrame)
    reduce(vcat, dfs) |> normalize_df! |> process_df!
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

    sort!(df, [:time, :z, :y, :x])
end