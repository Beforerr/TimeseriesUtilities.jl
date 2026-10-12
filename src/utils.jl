/ₜ(x, n) = x / n
/ₜ(x::P, n) where {P <: Dates.AbstractTime} = P(cld(Dates.value(x), n))

# https://github.com/JuliaStats/Statistics.jl/issues/203
@inline function _quantile_sorted(v, p::Real)
    n = length(v)
    n == 1 && @inbounds return v[1]
    h = (n - 1) * p
    i = floor(Int, h)
    f = h - i
    @inbounds return v[i + 1] * (1 - f) + v[min(i + 2, n)] * f
end

@inline function _median_sorted(v, lo::Int, hi::Int)
    n = hi - lo + 1
    m = lo + n ÷ 2
    return isodd(n) ? v[m] : middle(v[m - 1], v[m])
end

# `map` narrows a Union eltype to values seen, so result type depend on data; 
# use inferred eltype when it is concrete or a union of concrete types.
function _map(f::F, xs) where {F}
    T = Base.promote_op(f, eltype(xs))
    return _isconcreteunion(T) ? map!(f, similar(xs, T), xs) : map(f, xs)
end
_isconcreteunion(T) = T isa Union ? _isconcreteunion(T.a) && _isconcreteunion(T.b) : isconcretetype(T)

# static length, so slicing or rebuilding over the result infers when `d` is constant
other_dims(A, d) = ntuple(i -> i < d ? i : i + 1, Val(ndims(A) - 1))

# Assemble per-channel results (ordered as `vec` of the non-time dims) into `A`'s layout.
function _stack_channels(A, d, cols)
    out = similar(A, eltype(first(cols)), ntuple(i -> i == d ? length(first(cols)) : size(A, i), ndims(A)))
    out3 = _as3d(out, d)
    for (k, I) in enumerate(CartesianIndices((axes(out3, 1), axes(out3, 3))))
        out3[I[1], :, I[2]] = cols[k]
    end
    return out
end

# View `A` as (dims before `d`, `d`, dims after `d`): kernels then compile once per
# array type instead of once per (type, d).
_as3d(A, d) = reshape(A, prod(i -> size(A, i), 1:(d - 1); init = 1), size(A, d), :)

# Call `f(p, q)` over the non-time indices of an `_as3d` array. Time on the first/last dim
# leaves a length-1 loop whose overhead dominates few-channel data (e.g. 3-vectors);
@inline function _foreach_pq(f, A3)
    P, Q = size(A3, 1), size(A3, 3)
    if P == 1
        for q in 1:Q
            f(1, q)
        end
    elseif Q == 1
        for p in 1:P
            f(p, 1)
        end
    else
        for q in 1:Q, p in 1:P
            f(p, q)
        end
    end
    return
end

# https://github.com/joshday/SearchSortedNearest.jl
function searchsortednearest(a, x; by = identity, lt = isless, rev = false, distance = (a, b) -> abs(a - b))
    i = searchsortedfirst(a, x; by, lt, rev)
    if i == 1
    elseif i > length(a)
        i = length(a)
    elseif a[i] == x
    else
        i = lt(distance(by(a[i]), by(x)), distance(by(a[i - 1]), by(x))) ? i : i - 1
    end
    return i
end

_float_copy(A) = (data = unwrap(A); rebuild(A, copyto!(similar(data, float(eltype(data))), data)))

_check_holds_nan(data, f) = eltype(data) <: Union{AbstractFloat, Complex{<:AbstractFloat}} ||
    throw(ArgumentError("$(f)! writes NaN, which an array of $(eltype(data)) cannot hold; use $f"))
