using DataFrames,
    DataFramesMeta,
    CategoricalArrays
using Statistics

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