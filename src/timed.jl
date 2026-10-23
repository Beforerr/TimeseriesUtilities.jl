"""
    Timed(data, t; dim=nothing)
    Timed{D}(data, keys)

Plain array `data` with times `t` along dimension `dim` (optional for vectors); operations return a `Timed`.
`Timed{D}(data, keys)` takes keys for every axis, `D` being the time axis.
"""
struct Timed{D, A, K <: Tuple}
    data::A
    keys::K
    function Timed{D}(data::A, keys::K) where {D, A, K <: Tuple}
        D isa Int && 1 <= D <= ndims(data) || throw(ArgumentError("time dimension $D is out of range"))
        length(keys) == ndims(data) || throw(DimensionMismatch("$(length(keys)) keys for $(ndims(data)) dimensions"))
        isnothing(keys[D]) && throw(ArgumentError("time dimension $D has no times"))
        for (i, k) in enumerate(keys)
            isnothing(k) || _fitskey(k, data, i, D) ||
                throw(DimensionMismatch("$(summary(k)) keys do not fit dimension $i of data of size $(size(data))"))
        end
        return new{D, A, K}(data, keys)
    end
end

function Timed(data, t; dim = nothing)
    d = resolve_dim(data, dim)
    return Timed{Int(d)}(data, ntuple(i -> i == d ? t : nothing, ndims(data)))
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

_fitskey(k, data, i, D) = length(k) == size(data, i)
_fitskey(k::AbstractMatrix, data, i, D) = i != D && size(k) == (size(data, min(i, D)), size(data, max(i, D)))

rebuild(x::Timed{D}, data, dims::Tuple) where {D} =
    _timed(_findval(Val(D), map(_axis, dims)), data, map(d -> _newkey(x, d, dims), dims))

_newkey(_, d::Pair, _) = last(d)
_newkey(x::Timed{D}, ::Val{i}, dims) where {D, i} = _movekey(x.keys[i], Val(i < D), _findval(Val(i), dims), _findval(Val(D), dims))
_movekey(k, _, _, _) = k
# without the time axis as it was, a time-varying key no longer fits
_movekey(::AbstractMatrix, _, _, ::Nothing) = nothing
_movekey(k::AbstractMatrix, ::Val{before}, ::Val{j}, ::Val{p}) where {before, j, p} = before == (j < p) ? k : permutedims(k)

# same samples relabeled, so keys varying in time still fit, where `rebuild` assumes new samples and drops them
@inline rekey(x::Timed{D}, dim, keys) where {D} = Timed{D}(x.data, Base.setindex(x.keys, keys, dim))

Base.:(==)(x::Timed, y::Timed) = tdimnum(x) == tdimnum(y) && x.keys == y.keys && x.data == y.data
Base.isequal(x::Timed, y::Timed) = tdimnum(x) == tdimnum(y) && isequal(x.keys, y.keys) && isequal(x.data, y.data)
Base.hash(x::Timed, h::UInt) = hash(x.data, hash(x.keys, hash(tdimnum(x), hash(Timed, h))))

function Base.show(io::IO, x::Timed)
    print(io, "Timed(", summary(x.data), ", ")
    _show_span(io, times(x))
    print(io, "; dim = ", tdimnum(x), ")")
end

function Base.show(io::IO, m::MIME"text/plain", x::Timed)
    print(io, "Timed along dimension ", tdimnum(x), ", ")
    _show_span(io, times(x))
    println(io)
    show(io, m, x.data)
end

_show_span(io, t) = isempty(t) ? print(io, "no times") : print(io, length(t), " times ", first(t), " … ", last(t))

# Without a time axis the result is plain data.
_timed(::Nothing, data, _) = data
_timed(::Val{D}, data, keys) where {D} = Timed{D}(data, keys)

# `Val` position of `v` in `dims`, by dispatch: constant propagation through `===` gives up inside larger callers
_findval(v, dims::Tuple) = _findval(v, dims, Val(1))
_findval(::V, ::Tuple{V, Vararg}, i::Val) where {V} = i
_findval(v, dims::Tuple, ::Val{i}) where {i} = _findval(v, Base.tail(dims), Val(i + 1))
_findval(_, ::Tuple{}, ::Val) = nothing

Base.getindex(x::Timed, I...) = _index(getindex, x, I)
Base.view(x::Timed, I...) = _index(view, x, I)

# Scalar indices drop their dimension along with its keys.
function _index(f, x::Timed{D}, I) where {D}
    length(I) == ndims(x) || throw(ArgumentError("index a `Timed` with one index per dimension"))
    data = f(x.data, I...)
    I[D] isa Integer && return data
    # filtered by type so the keys infer: positions `j` would index `x.keys` at runtime
    keys = filter(!(k -> k isa _Dropped), map((k, i, j) -> i isa Integer ? _Dropped() : _slicekey(f, k, i, j, I, Val(D)), x.keys, I, ntuple(Val, Val(length(I)))))
    return Timed{D - count(i -> i isa Integer, I[1:(D - 1)])}(data, keys)
end

struct _Dropped end
_slicekey(f, k, i, _...) = f(k, i)
_slicekey(_, ::Nothing, _...) = nothing
_slicekey(f, k::AbstractMatrix, _, ::Val{j}, I, ::Val{D}) where {j, D} = f(k, I[min(j, D)], I[max(j, D)])
