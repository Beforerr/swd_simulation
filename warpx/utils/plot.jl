z_norm_lab = L"z~(d_i)"
time_norm_lab = L"T~(T_{ci})"

ids = [:z_norm, :time_norm]
labs = [z_norm_lab, time_norm_lab]

B_fields = ["Bx", "By", "Bz", "Bmag"]
E_fields = ["Ex", "Ey", "Ez", "Emag"]
j_field = ["jx", "jy", "jz", "jmag"]
temp_f_fields = ["T_parp", "T_perp"]
temp_norm_fields = ["T_parp_norm", "T_perp_norm"]

function plot_fields(df, fields; axis=NamedTuple(), figure=NamedTuple(), fig_options=(size=(800, 800),))

    temp_df = get_avg_fields(df, fields, ids=ids)
    plt = data(temp_df) * mapping(Pair.(ids, labs)..., :value, layout=:variable) * visual(Heatmap)

    draw(plt; figure=figure, axis=axis)
end

function plot_fields(df)

    pressure_fields = ["pressure_x", "pressure_y", "pressure_z"]
    pressure_f_fields = ["pressure_parp", "pressure_perp"]

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
    j_fields = names(df, r"j")

    df = select_time(df; window...)

    figure = (size=(600, 800),)

    axis = (ylabel="Current Density (A/m^2)",)
    plot_fields_time(df, j_fields; figure=figure, axis=axis)
    easy_save("j_time")

    axis = (ylabel="Magnetic Field (T)",)
    plot_fields_time(df, B_fields; figure=figure, axis=axis)
    easy_save("B_time")
end


function plot_field(df, field)
    temp_df = stack(df, field, ids)
    data(temp_df) * mapping(:z_norm => z_norm_lab, :value, color=:variable) * visual(Lines)
end

function plot_overview_ts(df;
    window=(; step=8),
)
    temp_df = select_time(df; window...)
    gdf = groupby(temp_df, :time_norm)

    fields = [B_fields, "rho_n_norm", j_field, temp_norm_fields, "anisotropy"]
    labels = ["B", "Normalized Density", "Current Density", "Temperature", "Anisotropy"]
    
    for (key, subdf) in pairs(gdf)
        fig = Figure(size=(1200, 1000),)
        axs = [Axis(fig[i, 1]) for i in 1:length(fields)]
    
        plts = map(plot_field $ subdf , fields)
        map(draw!, axs, plts)
        map(hidexdecorations!, axs)
    
        # add labels
        for (ax, label) in zip(axs, labels)
            ax.ylabel = label
        end
    
        axs[1].title = "Time: $(key.time_norm)"
        axs[end].xlabel = z_norm_lab
    
        easy_save("time_$(key.time_norm|>floor)", dir="figures/fields")
    end

end