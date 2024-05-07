z_norm_lab = L"z~(d_i)"
time_norm_lab = L"T~(T_{ci})"

ids = [:z_norm, :time_norm]
labs = [z_norm_lab, time_norm_lab]

function plot_fields(df, fields; axis=NamedTuple(), figure=NamedTuple(), fig_options=(size=(800, 800),))

    temp_df = get_avg_fields(df, fields, ids=ids)
    plt = data(temp_df) * mapping(Pair.(ids, labs)..., :value, layout=:variable) * visual(Heatmap)

    draw(plt; figure=figure, axis=axis)
end

function plot_fields(df)

    B_fields = names(df, r"^B(x|y|z|mag)$")
    E_fields = names(df, r"E")
    j_fields = names(df, r"j")
    pressure_fields = ["pressure_x", "pressure_y", "pressure_z"]
    pressure_f_fields = ["pressure_parp", "pressure_perp"]
    temp_f_fields = ["T_parp", "T_perp"]
    temp_norm_fields = ["T_parp_norm", "T_perp_norm"]


    plot_fields(df, B_fields)
    easy_save("B_field")
    plot_fields(df, E_fields)
    easy_save("E_field")
    plot_fields(df, j_fields)
    easy_save("j_field")

    plot_fields(df, "rho_n_norm")
    easy_save("rho_n_norm")

    plot_fields(df, temp_norm_fields)
    easy_save("temp_norm")

    plot_fields(df, "anisotropy")
    easy_save("anisotropy")
end