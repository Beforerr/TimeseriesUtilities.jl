# `(coords[i] - before, coords[i] + after)` for each index `i`.
struct PointWindows{T, C, B, A} <: AbstractVector{T}
    coords::C
    before::B
    after::A
end

function PointWindows(coords::C, before::B, after::A) where {C, B, A}
    T = Base.promote_op(_point_window, eltype(C), B, A)
    return PointWindows{T, C, B, A}(coords, before, after)
end

_point_window(t, before, after) = (t - before, t + after)

Base.size(pw::PointWindows) = (length(pw.coords),)
Base.getindex(pw::PointWindows, i::Int) = _point_window(pw.coords[i], pw.before, pw.after)


"""
    WindowedView(data, coords, windows; dim=1)
    WindowedView(data, coords, before, after; dim=1)

`AbstractVector` of views into `data`, one per time window. `wv[i]` returns the slice of `data` along dimension `dim` whose `coords` fall within the i-th window.

`windows` is any indexable source of `(t_start, t_stop)` pairs (e.g. `TimeRanges`). The `before`/`after` form builds a symmetric window around each point in `coords`. Requires `coords` to be sorted.
"""
struct WindowedView{dim, T, D, C, W} <: AbstractVector{T}
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

Base.size(wv::WindowedView) = (length(wv.windows),)

function Base.getindex(wv::WindowedView{dim}, i::Int) where {dim}
    t_start, t_stop = wv.windows[i]
    lo = searchsortedfirst(wv.coords, t_start)
    hi = searchsortedfirst(wv.coords, t_stop) - 1
    return _dim_view(wv.data, Val(dim), lo:hi)
end

"""
    twindows(x, every; period=every, start_by=:window, dim=nothing)

Windows of `x` as `ws, starts = twindows(x, every)`: `ws[i]` is a view of `x` over `[starts[i], starts[i] + period)`,
on the grid of [`tstat`](@ref); empty windows are kept as empty views.
"""
function twindows(x, every; period = every, start_by = :window, dim = nothing)
    d = dimnum(x, dim)
    t = axiskeys(x, d)
    windows = _tgrid_windows(t, every, period, start_by)
    return WindowedView{Int(d)}(x, t, windows), windows.coords
end

# Windows `[start, start + period)` with starts a range `every` apart.
function _tgrid_windows(t, every, period, start_by)
    isempty(t) && throw(ArgumentError("no samples to window"))
    issorted(t) || throw(ArgumentError("times must be sorted"))
    start_by in (:window, :datapoint) || throw(ArgumentError("start_by must be :window or :datapoint, got $(repr(start_by))"))
    t0 = start_by == :window ? _floor(first(t), every) : first(t)
    return PointWindows(t0:every:last(t), zero(period), period)
end

_floor(t, dt) = floor(t, dt)
_floor(i::Number, di) = i - mod(i, di)

"""
    smooth(x, window; dim=nothing, op=nanmean)

Moving `op` of `x` along `dim` over the windows `[c - before, c + after)` around each coordinate `c`;
a scalar `window` is split evenly into `before` and `after`.
Edge windows are truncated, so the output has the size of `x`.

Coordinates are the axis keys of `x`, or sample indices for a plain array.
"""
function smooth(x, window; dim = nothing, op = nanmean)
    d = dimnum(x, dim)
    windows = _point_windows(x, d, window)
    out = mapslices(unwrap(x); dims = d) do slice
        op.(WindowedView{1}(slice, windows.coords, windows))
    end
    return rebuild(x, out)
end

function _point_windows(x, d, window)
    coords = _haskeys(x) ? axiskeys(x, d) : axes(x, d)
    issorted(coords) || throw(ArgumentError("coordinates must be sorted"))
    return PointWindows(coords, _window_offsets(window, eltype(coords))...)
end

function _window_offsets(window::Tuple, _)
    @assert length(window) == 2
    return window[1], window[2]
end

_window_offsets(window, _) = window / 2, window / 2
_window_offsets(window::Real, T::Type{<:AbstractTime}) =
    throw(ArgumentError("a `$(typeof(window))` window does not apply to `$T` coordinates; pass a `Period`"))

# Split in the coordinates' resolution: on that grid `[c - w÷2, c + w - w÷2)` selects the same points as
# `[c - w/2, c + w/2)`, whereas a finer half-width would be silently rounded by e.g. `DateTime ± Nanosecond`.
function _window_offsets(window::Period, T)
    P = Base.promote_op(-, T, T)
    w = P <: Period && isconcretetype(P) ? convert(P, window) : Nanosecond(window)
    before = w ÷ 2
    return before, w - before
end
