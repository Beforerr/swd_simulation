z_norm_lab = L"z~(d_i)"
time_norm_lab = L"T~(T_{ci})"

ids = [:z_norm, :time_norm]
labs = [z_norm_lab, time_norm_lab]

B_fields = ["Bx", "By", "Bz", "Bmag"]
E_fields = ["Ex", "Ey", "Ez", "Emag"]
j_fields = ["jx", "jy", "jz", "jmag"]
temp_f_fields = ["T_parp", "T_perp"]
temp_norm_fields = ["T_parp_norm", "T_perp_norm"]
pressure_fields = ["pressure_x", "pressure_y", "pressure_z"]
pressure_f_fields = ["pressure_parp", "pressure_perp"]

function plot_fields(df, fields; axis=NamedTuple(), figure=NamedTuple(), fig_options=(size=(800, 800),))

    temp_df = get_avg_fields(df, fields, ids=ids)
    plt = data(temp_df) * mapping(Pair.(ids, labs)..., :value, layout=:variable) * visual(Heatmap)

    draw(plt; figure=figure, axis=axis)
end

function plot_fields(df)

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

    plot_fields(df, "Λ_temp")
    easy_save("Λ_temp")
end

function plot_fields_time(
    df, fields;
    axis=NamedTuple(), figure=NamedTuple(), legend=(framevisible=false,)
)
    temp_df = get_avg_fields(df, fields, ids=ids)
    temp_df.time_norm = CategoricalArray(temp_df.time_norm .|> floor)

    plt = data(temp_df) * mapping(:z_norm => z_norm_lab, :value, color=:variable, row=:time_norm) * visual(Lines)
    draw(plt; axis=axis, figure=figure, legend=legend)
end

function plot_fields_time(df; window=(; step=16))
    B_fields = names(df, r"^B(x|y|mag)$")

    df = select_time(df; window...)

    figure = (size=(600, 800),)

    axis = (ylabel="Current Density (A/m^2)",)
    plot_fields_time(df, j_fields; figure=figure, axis=axis)
    easy_save("j_time")

    axis = (ylabel="Magnetic Field (T)",)
    plot_fields_time(df, B_fields; figure=figure, axis=axis)
    easy_save("B_time")
end

function plot_field_1d(df; fields, ax, ylabel)
    lns = [lines!(ax, df.z_norm, df[:, field] .|> ustrip, label=field) for field in fields]

    axislegend(ax; framevisible=false)
    ax.ylabel = ylabel
    return lns
end

plot_B = plot_field_1d$(; fields=B_fields, ylabel="Magnetic Field (T)")
plot_j = plot_field_1d$(; fields=j_fields, ylabel="Current Density (A/m^2)")
plot_density = plot_field_1d$(; fields=["rho_n_norm"], ylabel="Normalized Density")
plot_anisotropy = plot_field_1d$(; fields=["Λ"], ylabel="Anisotropy")
plot_velocity = plot_field_1d$(; fields=["velocity_y", "vA_y", "vA_p_y"], ylabel="Velocity (m/s)")
plot_temp = plot_field_1d$(; fields=["T_z_norm"; temp_norm_fields], ylabel="Normalized Temperature")


function plot_overview_ts(df;
    window=(; step=8),
    plot_funcs=(plot_B, plot_j, plot_temp, plot_anisotropy, plot_velocity, plot_density)
)
    temp_df = select_time(df; window...)
    gdf = groupby(temp_df, :time_norm)

    for (key, subdf) in pairs(gdf)
        fig = Figure(size=(1200, 1000),)
        axs = [Axis(fig[i, 1]) for i in 1:length(plot_funcs)]
        fgs = [plot_func(subdf; ax=ax) for (ax, plot_func) in zip(axs, plot_funcs)]

        map(hidexdecorations!, axs[1:end-1])
        axs[1].title = "Time: $(key.time_norm|>round)"
        axs[end].xlabel = z_norm_lab

        easy_save("time_$(key.time_norm|>round)", dir="figures/fields")
    end

end