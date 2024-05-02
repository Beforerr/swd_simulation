using DataFrames,
    DataFramesMeta,
    CategoricalArrays

function get_avg_fields(df, fields; ids=ids)
    @chain df begin
        groupby(ids)
        combine(fields .=> mean, renamecols=false)
        stack(fields, ids)
    end
end

function select_time(df; start=1, stop=missing, step=1)
    gdf = groupby(df, :time)
    if ismissing(stop)
        stop = lastindex(gdf)
    end
    combine(gdf[start:step:stop], names(df))
end

function plot_fields_time(df, fields; ids=ids, window=missing, norm=false)

    if !ismissing(window)
        df = select_time(df; window...)
    end

    temp_df = get_avg_fields(df, fields, ids=ids)
    temp_df.time_norm = CategoricalArray(temp_df.time_norm .|> floor)

    if norm
        temp_df = @by temp_df :variable begin
            :value = :value / maximum(abs.(:value))
            :time_norm
            :z_norm
        end
    end

    # height = 200 * length(unique(temp_df.time_norm))
    data(temp_df) * mapping(:z_norm => z_norm_lab, :value, color=:variable, row=:time_norm) * visual(Lines)
end
