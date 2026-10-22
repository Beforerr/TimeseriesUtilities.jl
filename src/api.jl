# Backend interface. A container backend (an extension, or `Timed` for plain arrays) implements
# `tdimnum`, `dimnum` for its dimension-name types, `axiskeys`, `unwrap` and `rebuild`,
# plus indexing (`getindex`/`view` with one index per dimension) that keeps the keys of the kept axes.
# Backends whose lookups record sort order also implement `sorted_axis`.

"""
    tdimnum(x)

Index of the time dimension of `x`, or `nothing` if `x` has none; plain arrays have none.
"""
tdimnum(x) = nothing

"""
    dimnum(x, dim)

Index of dimension `dim` in `x`. Containers with named dimensions add methods for their name types.
"""
dimnum(_, dim::Integer) = dim

"""
    resolve_dim(x, dim = nothing)

Index of the dimension an operation's `dim` keyword selects in `x`.
"""
resolve_dim(x, dim) = dimnum(x, dim)
function resolve_dim(x, ::Nothing = nothing)
    d = tdimnum(x)
    isnothing(d) || return d
    # Containers declare their time dimension; a keyed vector without one is not a time series.
    x isa AbstractVector && !_iscontainer(x) && return 1
    throw(ArgumentError("$(typeof(x)) has no time dimension; pass `dim`"))
end

"""
    axiskeys(x, dim::Integer)

Coordinate values along dimension `dim` of `x`, or `nothing` without them (plain arrays, axes without coordinates).
Coordinates that vary in time are a matrix along `dim` and the time dimension, in the order of the dimensions of `x`.
"""
axiskeys(_, _) = nothing

"""
    unwrap(x)

Underlying array of container `x` (sharing memory), or `x` itself when it is not a container.
"""
unwrap(x) = x

# compared by type so it folds at compile time
_iscontainer(x) = typeof(unwrap(x)) !== typeof(x)

"""
    rebuild(x, data, dims::Tuple)
    rebuild(x, data, dim::Integer, keys)
    rebuild(x, data)

Container like `x` holding `data`, with axes `dims`: each `Val(i)` (dimension `i` of `x` as it is,
with whatever the container attaches to it), `Val(i) => keys` (dimension `i` with new samples),
`name::Symbol => keys` (new named dimension) or `nothing => keys` (new anonymous dimension).
`keys` is `nothing` for an axis without coordinates;
backends that cannot represent that supply positions.

Backends implement the 3-argument form, and the shorthands only where that is cheaper.
Plain arrays rebuild to `data`; containers without a method throw rather than silently drop their keys.
"""
rebuild(x, data) = rebuild(x, data, ntuple(Val, Val(ndims(x))))
@inline rebuild(x, data, dim, keys) = rebuild(x, data, Base.setindex(ntuple(Val, Val(ndims(x))), Val(dim) => keys, dim))
rebuild(x, data, _) = _iscontainer(x) ? throw(ArgumentError("$(typeof(x)) has no `rebuild` method")) : data
_axis(d::Val) = d
_axis(d::Pair) = first(d)
_axiskeys(x, ::Val{i}) where {i} = axiskeys(x, i)
_axiskeys(_, d::Pair) = last(d)

"""
    rekey(x, dim::Integer, keys)

`x` with the keys of dimension `dim` replaced.
"""
@inline rekey(x, dim, keys) =
    _iscontainer(x) ? rebuild(x, unwrap(x), dim, keys) : throw(ArgumentError("$(typeof(x)) has no axis keys to replace"))

# Hook for backends whose lookups record sort order (DimensionalData).
sorted_axis(sorted, dim; rev = false) = sorted

"""
    times(x, dim=nothing)

Time coordinate of `x`: the keys of its time dimension (or of `dim`); throws without them.
"""
@inline function times(x, dim = nothing)
    d = resolve_dim(x, dim)
    return @something axiskeys(x, d) throw(ArgumentError("dimension $d of $(typeof(x)) has no times; wrap it with `Timed(x, t)`"))
end

@inline _keys_or_indices(x, d) = @something axiskeys(x, d) axes(x, d)

# Apply `f(data, t, d)` to the underlying data of `A` and the keys `t` of dimension `d`,
# and rebuild `A` from the returned `(data, t)`.
@inline function _ontime(f::F, A, dim) where {F}
    d = resolve_dim(A, dim)
    out, t = f(unwrap(A), times(A, d), d)
    return rebuild(A, out, d, t)
end

_valof(::Val{D}) where {D} = D

# Rebuild a reduction over dimension `d` of `x` with the remaining axes.
@inline _dropdim(x, data, d) = _iscontainer(x) ? _dropdim(x, data, Val(d)) : data
_dropdim(x, data, ::Val{d}) where {d} = rebuild(x, data, map(Val, other_dims(x, d)))

_seconds(dt::Period) = dt / Second(1)
_seconds(dt) = dt
