"""
    eachsample(A, [ts]; dim=nothing)

Iterate `(x, t)`: each slice `x` of `A` along time with its time `t`.
"""
eachsample(A; dim = nothing) = eachsample(A, axiskeys(A, dimnum(A, dim)); dim)
eachsample(A, ts; dim = nothing) = (d = dimnum(A, dim); zip(eachslice(unwrap(A); dims = d), _sampletimes(A, ts, d)))

@inline function _sampletimes(A, ts, d)
    bt = Base.broadcastable(ts)
    N = length(axes(bt))  # `ndims` is undefined for tuples
    N == 0 && return Iterators.repeated(bt[])
    N == 1 || throw(ArgumentError("ts must be 0- or 1-dimensional, got $N dimensions"))
    length(bt) == size(A, d) ||
        throw(DimensionMismatch("length(ts) = $(length(bt)) must match size(A, $d) = $(size(A, d))"))
    return bt
end

"""
    mapsamples(f, A, [ts]; dim=nothing)

Apply `f(x, t)` per sample and stack results along time.

Results the size of a sample keep the other axes of `A`; other results get anonymous axes.
"""
@inline mapsamples(f, A; dim = nothing) = mapsamples(f, A, axiskeys(A, dimnum(A, dim)); dim)
@inline mapsamples(f, A, ts; dim = nothing) =
    (d = dimnum(A, dim); _mapsamples(f, A, _sampletimes(A, ts, d), d, Val(ndims(A))))

# `eachslice` and the rebuild need a static `D`: `Val(d)` with a runtime `d` infers `Any`, testing `d == N`
# for each `N` union-splits, and a constant `d` (default `dim`, wrapper dimension types) folds to one branch.
@inline _mapsamples(f, A, ts, d, D::Val{N}) where {N} =
    d == N ? _rebuild_samples(A, _stack_samples(f, unwrap(A), ts, D), D) : _mapsamples(f, A, ts, d, Val(N - 1))
_mapsamples(_, _, _, d, ::Val{0}) = throw(ArgumentError("dim = $d is out of range"))

function _stack_samples(f, A, ts, ::Val{D}) where {D}
    xs = eachslice(A; dims = D)
    isempty(xs) && return _empty_samples(A, Base.promote_op(f, eltype(xs), eltype(ts)), D)
    samples = zip(xs, ts)
    (x, t), state = iterate(samples)
    y = f(x, t)
    P = _slotdim(y, D)
    return _fill_samples!(_similar_samples(A, y, length(xs), P), f, y, samples, state, 1, P)
end

# `nothing`: each result is an element of `out`; `Val(p)`: each result fills slice `k` along axis `p`
_slotdim(_, _) = nothing
_slotdim(y::AbstractArray, d) = Val(_timepos(d, ndims(y)))

# position of the time axis in the output: `d`, clamped when result has fewer axes than a sample
_timepos(d, nresult) = min(d, nresult + 1)

_similar_samples(A, y, n, ::Nothing) = similar(A, typeof(y), n)
function _similar_samples(A, y, n, ::Val{P}) where {P}
    sz = size(y)
    return similar(A, eltype(y), ntuple(i -> i < P ? sz[i] : i == P ? n : sz[i - 1], Val(ndims(y) + 1)))
end

# no sample: result type comes from inference and per-sample size from `A`
_empty_samples(A, T, _) = similar(A, T, 0)
function _empty_samples(A, ::Type{<:AbstractArray{S, M}}, d) where {S, M}
    M == ndims(A) - 1 && return similar(A, S, Base.setindex(size(A), 0, d))
    return similar(A, S, ntuple(Returns(0), M + 1))
end

# When its element type does not fit, widen `out` and continue from `k`
# (same strategy as `map`; for a type-stable `f` the check folds away).
function _fill_samples!(out, f, y, samples, state, k, P)
    while true
        T = isnothing(P) ? typeof(y) : eltype(y)
        if !(T <: eltype(out))
            new = similar(out, Base.promote_typejoin(eltype(out), T))
            for j in 1:(k - 1)
                _setsample!(new, j, _getsample(out, j, P), P)
            end
            return _fill_samples!(new, f, y, samples, state, k, P)
        end
        _setsample!(out, k, y, P)
        next = iterate(samples, state)
        isnothing(next) && return out
        (x, t), state = next
        y, k = f(x, t), k + 1
    end
end

_getsample(out, k, ::Nothing) = out[k]
_getsample(out, k, ::Val{P}) where {P} = view(out, ntuple(Returns(:), P - 1)..., k, ntuple(Returns(:), ndims(out) - P)...)

_setsample!(out, k, y, ::Nothing) = (out[k] = y)

function _setsample!(out, k, y, ::Val{P}) where {P}
    sz = size(out)
    slot = ntuple(i -> sz[i < P ? i : i + 1], Val(ndims(out) - 1))
    size(y) == slot || throw(DimensionMismatch("sample $k result has size $(size(y)); the first had $slot"))
    # explicit loops so static `y` can unroll here for vectorization
    if P == ndims(out)  # time last: slots are contiguous
        off = (k - 1) * length(y)
        @inbounds for (j, v) in enumerate(y)
            out[off + j] = v
        end
    else
        @inbounds for (I, v) in zip(CartesianIndices(slot), y)
            out[Tuple(I)[1:(P - 1)]..., k, Tuple(I)[P:end]...] = v
        end
    end
    return out
end

function _rebuild_samples(A, data, ::Val{D}) where {D}
    unwrap(A) === A && return data
    size(data) == size(A) && return rebuild(A, data)
    p = _timepos(D, ndims(data) - 1)
    newdims = Base.setindex(ntuple(Returns(nothing), Val(ndims(data))), Val(D), p)
    return rebuild(A, data, newdims, Base.setindex(axes(data), axiskeys(A, D), p))
end
