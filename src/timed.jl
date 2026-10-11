"""
    Timed(data, t; dim=nothing)
    Timed{D}(data, keys)

Plain array `data` with times `t` along dimension `dim` (optional for vectors); operations return a `Timed`.
`Timed{D}(data, keys)` takes keys for every axis, `D` being the time axis.

Supports this package's operations and indexing with one index per dimension. It is not an `AbstractArray`,
so other array code (broadcasting, `sum`, plotting) throws rather than silently dropping the times;
take it apart first with `data, t = x`.
"""
struct Timed{D, A, K <: Tuple}
    data::A
    keys::K
    function Timed{D}(data::A, keys::K) where {D, A, K <: Tuple}
        D isa Int && 1 <= D <= ndims(data) || throw(ArgumentError("time dimension $D is out of range"))
        length(keys) == ndims(data) || throw(DimensionMismatch("$(length(keys)) keys for $(ndims(data)) dimensions"))
        for (i, k) in enumerate(keys)
            length(k) == size(data, i) ||
                throw(DimensionMismatch("$(length(k)) keys for dimension $i of size $(size(data, i))"))
        end
        return new{D, A, K}(data, keys)
    end
end

function Timed(data, t; dim = nothing)
    d = dimnum(data, dim)
    return Timed{Int(d)}(data, ntuple(i -> i == d ? t : axes(data, i), ndims(data)))
end

# `data, t = x` lowers to `indexed_iterate`; defining `iterate` instead would make `Timed` a collection.
Base.indexed_iterate(x::Timed, i::Int, state = 1) =
    i == 1 ? (x.data, 2) : i == 2 ? (times(x), 3) : throw(BoundsError(x, i))

Base.parent(x::Timed) = x.data
Base.size(x::Timed, d...) = size(x.data, d...)
Base.axes(x::Timed, d...) = axes(x.data, d...)
Base.ndims(x::Timed) = ndims(x.data)
Base.eltype(::Type{<:Timed{D, A}}) where {D, A} = eltype(A)

unwrap(x::Timed) = x.data
tdimnum(::Timed{D}) where {D} = D
axiskeys(x::Timed, d::Integer) = x.keys[d]

rebuild(x::Timed{D}, data) where {D} = Timed{D}(data, x.keys)
rebuild(x::Timed{D}, data, d::Integer, keys) where {D} = Timed{D}(data, Base.setindex(x.keys, keys, d))
rebuild(x::Timed{D}, data, dims::Tuple, keys::Tuple) where {D} = _timed(_findval(Val(D), dims), data, keys)
rekey(x::Timed, d, keys) = rebuild(x, x.data, d, keys)

Base.:(==)(x::Timed, y::Timed) = tdimnum(x) == tdimnum(y) && x.keys == y.keys && x.data == y.data

# Without a time axis the result is plain data.
_timed(::Nothing, data, _) = data
_timed(D::Int, data, keys) = Timed{D}(data, keys)

# position of `v` in `dims`, by recursion so it folds to a constant for static `dims`
@inline _findval(v, dims::Tuple, i = 1) = first(dims) === v ? i : _findval(v, Base.tail(dims), i + 1)
_findval(v, ::Tuple{}, i = 1) = nothing

Base.getindex(x::Timed, I...) = _index(getindex, x, I)
Base.view(x::Timed, I...) = _index(view, x, I)

# Scalar indices drop their dimension along with its keys.
function _index(f, x::Timed{D}, I) where {D}
    length(I) == ndims(x) || throw(ArgumentError("index a `Timed` with one index per dimension"))
    data = f(x.data, I...)
    kept = filter(i -> !(I[i] isa Integer), ntuple(identity, length(I)))
    return _timed(_findval(D, kept), data, map(i -> f(x.keys[i], I[i]), kept))
end
