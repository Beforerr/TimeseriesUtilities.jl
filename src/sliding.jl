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
