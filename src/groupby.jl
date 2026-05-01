# https://dataframes.juliadata.org/stable/man/split_apply_combine/
# https://github.com/JuliaData/SplitApplyCombine.jl

# https://docs.pola.rs/api/python/stable/reference/dataframe/api/polars.DataFrame.group_by_dynamic.html

_floor(t, dt) = floor(t, dt)
_floor(i::Number, di) = i - mod(i, di)

"""
    GroupByDynamic(times, every, period=every, start_by=:window)

Iterator over `(index_range, window_start)` pairs for sliding windows on sorted `times`.
"""
struct GroupByDynamic{T, E, P}
    times::T
    every::E
    period::P
    start_by::Symbol
end

GroupByDynamic(times, every, period=every, start_by=:window) =
    GroupByDynamic(times, every, period, start_by)

Base.IteratorSize(::Type{<:GroupByDynamic}) = Base.SizeUnknown()
Base.eltype(::Type{<:GroupByDynamic{T}}) where {T} = Tuple{UnitRange{Int}, eltype(T)}

function Base.iterate(iter::GroupByDynamic, current_start = _initial_start(iter))
    (; times, every, period) = iter
    isempty(times) && return nothing
    max_t = last(times)
    current_start > max_t && return nothing
    while current_start <= max_t
        window_end = current_start + period
        start_idx = searchsortedfirst(times, current_start)
        end_idx = searchsortedfirst(times, window_end) - 1
        next_start = current_start + every
        start_idx <= end_idx && return (start_idx:end_idx, current_start), next_start
        current_start = next_start
    end
    return nothing
end

_initial_start(iter::GroupByDynamic) =
    ifelse(iter.start_by == :window, _floor(first(iter.times), iter.every), first(iter.times))

"""
    groupby_dynamic(x, every, period=every, start_by=:window)

Eagerly collect `(group_idx, starts)` vectors for sorted `x`.
"""
function groupby_dynamic(x, every, args...)
    iter = GroupByDynamic(x, every, args...)
    n = Base.min(floor(Int, (last(x) - _initial_start(iter)) / every) + 1, length(x))
    group_idx = Vector{UnitRange{Int}}(undef, n)
    starts = Vector{eltype(x)}(undef, n)
    i = 0
    for (idx, t) in iter
        i += 1
        group_idx[i] = idx
        starts[i] = t
    end
    resize!(group_idx, i)
    resize!(starts, i)
    return group_idx, starts
end


"""
    tgroupby(x, every, period = every, start_by = :window)

Group `x` into windows based on `every` and `period`.
"""
function tgroupby(x, args...; dim=nothing, kwargs...)
    d = dimnum(x, dim)
    times = axiskeys(x, d)
    group_idx, = groupby_dynamic(times, args...; kwargs...)
    return map(group_idx) do idx
        selectdim(x, d, idx)
    end
end
