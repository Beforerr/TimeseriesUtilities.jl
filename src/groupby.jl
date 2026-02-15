# https://dataframes.juliadata.org/stable/man/split_apply_combine/
# https://github.com/JuliaData/SplitApplyCombine.jl

# https://docs.pola.rs/api/python/stable/reference/dataframe/api/polars.DataFrame.group_by_dynamic.html

_floor(t, dt) = floor(t, dt)
_floor(i::Number, di) = i - mod(i, di)

"""
Group `x` into windows based on `every` and `period`.
"""
function groupby_dynamic(x, every, period = every, start_by = :window)
    min, max = timerange(x)
    current_start = ifelse(start_by == :window, _floor(min, every), min)
    n = Base.min(floor(Int, (max - current_start) / every) + 1, length(x))
    group_idx = Vector{UnitRange{Int}}(undef, n)
    starts = Vector{eltype(x)}(undef, n)
    i = 0
    while current_start <= max
        window_end = current_start + period
        # Find indices of rows that fall in the current window using searchsorted for better performance
        start_idx = searchsortedfirst(x, current_start)
        end_idx = searchsortedfirst(x, window_end) - 1
        if start_idx <= end_idx
            i += 1
            starts[i] = current_start
            group_idx[i] = start_idx:end_idx
        end
        current_start += every
    end
    resize!(group_idx, i)
    resize!(starts, i)
    return group_idx, starts
end


"""
    tgroupby(data, times, every, args...; dim=1, kwargs...)

Group `data` into windows based on `every` and `period` using `times` for grouping.
Returns a vector of array slices.
"""
function tgroupby(data::AbstractArray, ts::AbstractVector, args...; dim = 1, kwargs...)
    group_idx, = groupby_dynamic(ts, args...; kwargs...)
    return map(group_idx) do idx
        selectdim(data, dim, idx)
    end
end
