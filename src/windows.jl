abstract type AbstractWindowView{T} <: AbstractVector{T} end

Base.size(wv::AbstractWindowView) = (length(wv.windows),)

# Window spec yielding `(coords[i] - before, coords[i] + after)` for each index `i`.
struct PointWindows{C, B, A}
    coords::C
    before::B
    after::A
end

Base.length(pw::PointWindows) = length(pw.coords)
Base.eltype(::Type{PointWindows{C, B, A}}) where {C, B, A} = NTuple{2, eltype(C)}

function Base.getindex(pw::PointWindows, i::Int)
    t = pw.coords[i]
    return t - pw.before, t + pw.after
end


"""
    WindowedView(data, coords, windows; dim=1)
    WindowedView(data, coords, before, after; dim=1)

`AbstractVector` of views into `data`, one per time window. `wv[i]` returns the slice of `data` along dimension `dim` whose `coords` fall within the i-th window.

`windows` is any indexable source of `(t_start, t_stop)` pairs (e.g. `TimeRanges`). The `before`/`after` form builds a symmetric window around each point in `coords`. Requires `coords` to be sorted.
"""
struct WindowedView{dim, T, D, C, W} <: AbstractWindowView{T}
    data::D
    coords::C
    windows::W
end

_dim_view(data, ::Val{dim}, i) where {dim} =
    view(data, ntuple(k -> k == dim ? i : Colon(), ndims(data))...)
_dim_view(data::AbstractVector, ::Val{1}, i) = view(data, i)

function WindowedView{dim}(data::D, coords::C, windows::W) where {dim, D, C, W}
    T = Base.promote_op(_dim_view, D, Val{dim}, UnitRange{Int})
    return WindowedView{dim, T, D, C, W}(data, coords, windows)
end

function Base.getindex(wv::WindowedView{dim}, i::Int) where {dim}
    t_start, t_stop = wv.windows[i]
    lo = searchsortedfirst(wv.coords, t_start)
    hi = searchsortedfirst(wv.coords, t_stop) - 1
    return _dim_view(wv.data, Val(dim), lo:hi)
end

"""
    smooth(data, times, window; dim=nothing)
    smooth(data, window; dim=nothing)

Smooths a time series by computing a moving average over a sliding window.
Edge windows are truncated, so the output has the same size as the input.

Each window covers the half-open interval `[coord - before, coord + after)`.
A scalar `window` is interpreted as a coordinate span along the smoothed axis.

# Arguments
- `dim`: Dimension to smooth along; defaults to the time dimension
- `op=nanmean`: Function used to aggregate each window
"""
@inline function smooth(data, coords, window; dim = nothing, op = nanmean)
    dim = dimnum(data, dim)
    length(coords) == size(data, dim) || throw(DimensionMismatch("length(coords) must match size(data, dim)"))
    issorted(coords) || throw(ArgumentError("coords must be sorted"))
    before, after = _window_offsets(window, eltype(coords))
    windows = PointWindows(coords, before, after)
    return mapslices(data; dims = dim) do slice
        op.(WindowedView{1}(slice, coords, windows))
    end
end

function smooth(data, window; dim = nothing, kw...)
    d = dimnum(data, dim)
    return smooth(data, axes(data, d), window; dim = d, kw...)
end


function _window_offsets(window::Tuple, _)
    @assert length(window) == 2
    return window[1], window[2]
end

_window_offsets(window, _) = window / 2, window / 2

# Split in the coordinates' resolution: on that grid `[c - w÷2, c + w - w÷2)` selects the same points as
# `[c - w/2, c + w/2)`, whereas a finer half-width would be silently rounded by e.g. `DateTime ± Nanosecond`.
function _window_offsets(window::Period, T)
    P = Base.promote_op(-, T, T)
    w = P <: Period && isconcretetype(P) ? convert(P, window) : Nanosecond(window)
    before = w ÷ 2
    return before, w - before
end

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
