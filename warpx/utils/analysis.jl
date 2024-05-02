using DataFrames,
    DataFramesMeta,
    CategoricalArrays

function normalize_df!(df)
    @transform!(df,
        :time_norm = :time ./ meta["t_ci"],
        :z_norm = :z ./ meta["d_i"],
    )
end

function process_df!(df)
    normalize_df!(df)

    B_fields = names(df, r"B")
    E_field = names(df, r"E")
    j_field = names(df, r"j")

    transform!(df,
        B_fields => ByRow(norm ∘ vcat) => :Bmag,
        E_field => ByRow(norm ∘ vcat) => :Emag,
        j_field => ByRow(norm ∘ vcat) => :jmag,
    )
end


function plot_fields_time(df, fields; ids=ids, step=1, norm=false)

    # select subset of the data with step
    # TODO: implement this in a more efficient way
    gdf = groupby(df, :time)[begin:step:end]
    temp_df = combine(gdf, names(df))

    temp_df = get_avg_fields(temp_df, fields, ids=ids)
    temp_df.time_norm = CategoricalArray(temp_df.time_norm .|> floor)

    if norm
        temp_df = @by temp_df :variable begin
            :value = :value / maximum(abs.(:value))
            :time_norm
            :z_norm
        end
    end

    height = 100 * length(unique(temp_df.time_norm))
    fig_options = (size=(1200, height),)

    plt = data(temp_df) * mapping(:z_norm => z_norm_lab, :value, color=:variable, row=:time_norm) * visual(Lines)
    draw(plt; figure=fig_options)
end
