# Backend interface. A container backend (an extension, or `Timed` for plain arrays) implements
# `tdimnum`, `dimnum` for its dimension-name types, `axiskeys`, `unwrap` and `rebuild`,
# plus indexing (`getindex`/`view` with one index per dimension) that keeps the keys of the kept axes.
# Backends whose lookups record sort order also implement `sorted_axis`.

"""
    tdimnum(x)

Index of the time dimension of `x`, or `nothing` if `x` has none.
A vector's time dimension is assumed to be its only one, multidimensional arrays have no default time dimension.
"""
tdimnum(x) = nothing
tdimnum(x::AbstractArray) = ndims(x) == 1 ? 1 : nothing

"""
    dimnum(x, dim)

Index of dimension `dim` in `x`; `nothing` selects the time dimension ([`tdimnum`](@ref)).

Containers with named dimensions add methods for their name types.
"""
dimnum(_, dim::Integer) = dim
dimnum(x, ::Nothing) =
    @something tdimnum(x) throw(ArgumentError("$(typeof(x)) has no time dimension; pass `dim`"))

"""
    axiskeys(x, dim::Integer)

Coordinate values along dimension `dim` of a keyed container.
"""
axiskeys(x, _) = throw(ArgumentError("$(typeof(x)) has no axis keys; wrap it with `Timed(x, t)`"))

"""
    unwrap(x)

Underlying array of container `x` (sharing memory), or `x` itself when it is not a container.
"""
unwrap(x) = x

_haskeys(x) = unwrap(x) !== x

"""
    rebuild(x, data)
    rebuild(x, data, dim::Integer, keys)
    rebuild(x, data, dims::Tuple, keys::Tuple)

Container like `x` holding `data`: with the axes of `x`; with the keys of dimension `dim` replaced;
or with new axes `dims`, each a `Val(i)` (dimension `i` of `x`, static so the result is inferable),
a `Symbol` (new named dimension) or `nothing` (new anonymous dimension), with coordinates `keys`.

Plain arrays rebuild to `data`; containers without a method throw rather than silently drop their keys.
"""
rebuild(x, data) = _norebuild(x, data)
rebuild(x, data, _, _) = _norebuild(x, data)
_norebuild(x, data) = _haskeys(x) ? throw(ArgumentError("$(typeof(x)) has no `rebuild` method")) : data

"""
    rekey(x, dim::Integer, keys)

`x` with the keys of dimension `dim` replaced.
"""
@inline rekey(x, dim, keys) =
    _haskeys(x) ? rebuild(x, unwrap(x), dim, keys) : throw(ArgumentError("$(typeof(x)) has no axis keys to replace"))

# Hook for backends whose lookups record sort order (DimensionalData).
sorted_axis(sorted, dim; rev = false) = sorted

"""
    times(x, dim=nothing)

Time coordinate of `x`: the keys of its time dimension (or of `dim`).
"""
times(x, dim = nothing) = axiskeys(x, dimnum(x, dim))

# Apply `f(data, t, d)` to the underlying data of `A` and the keys `t` of dimension `d`,
# and rebuild `A` from the returned `(data, t)`.
@inline function _ontime(f::F, A, dim) where {F}
    d = dimnum(A, dim)
    out, t = f(unwrap(A), axiskeys(A, d), d)
    return rebuild(A, out, d, t)
end

_valof(::Val{D}) where {D} = D

# Rebuild a reduction over dimension `d` of `x` with the remaining axes.
@inline _dropdim(x, data, d) = _haskeys(x) ? _dropdim(x, data, Val(d)) : data
function _dropdim(x, data, ::Val{d}) where {d}
    others = map(Val, other_dims(x, d))
    return rebuild(x, data, others, map(v -> axiskeys(x, _valof(v)), others))
end

_seconds(dt::Period) = dt / Second(1)
_seconds(dt) = dt
