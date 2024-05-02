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

    B_fields = names(df, r"B")
    E_field = names(df, r"E")
    j_field = names(df, r"j")

    normalize_df!(df)
    transform!(df,
        B_fields => ByRow(norm ∘ vcat) => :Bmag,
        E_field => ByRow(norm ∘ vcat) => :Emag,
        j_field => ByRow(norm ∘ vcat) => :jmag,
    )
end