"""
    tsubtract(x, op=nanmedian; dim=nothing)

Subtract a statistic `op` along time dimension (or `dim`) from `x`.
"""
function tsubtract(x, op = nanmedian; dim = nothing)
    dims = dimnum(x, dim)
    return x .- op(parent(x); dims)
end

_deriv_tfunc(T1::Type, T2::Type) = identity
_deriv_tfunc(::Type{T1}, ::Type{T2}) where {T1 <: Real, T2 <: Dates.TimeType} = _seconds

_deriv_tfunc(A, t) = _deriv_tfunc(eltype(A), eltype(t))

function _tderiv(A, times; dim = 1)
    # return diff(A; dims) ./ diff(times) # this allocates and is slow
    Base.require_one_based_indexing(A)
    N = ndims(A)
    1 <= dim <= N || throw(ArgumentError("dimension $dim out of range (1:$N)"))

    r = Base.axes(A)
    r0 = ntuple(i -> i == dim ? UnitRange(1, last(r[i]) - 1) : UnitRange(r[i]), N)
    r1 = ntuple(i -> i == dim ? UnitRange(2, last(r[i])) : UnitRange(r[i]), N)
    rt0 = r0[dim]
    rt1 = r1[dim]
    tfunc = _deriv_tfunc(A, times)
    return (view(A, r1...) .- view(A, r0...)) ./ tfunc.(view(times, rt1) .- view(times, rt0))
end

"""
    tderiv(A, times; dim=1)
    tderiv(A; dim=nothing)

Compute the time derivative of `A`. Set `lazy=true` for lazy evaluation.

See also: [deriv_data - PySPEDAS](https://pyspedas.readthedocs.io/en/latest/_modules/pyspedas/analysis/deriv_data.html)
"""
@inline function tderiv(A, times = nothing; lazy = false, dim = nothing, kw...)
    dim = dimnum(A, dim)
    f = lazy ? DiffQ : _tderiv
    times = @something times axiskeys(A, dim)
    out = f(unwrap(A), times; dim)
    return rebuild(A, out, dim, @view times[1:(end - 1)])
end

"""
    proj(a, b)

Vector projection of a vector `a` on (or onto) a nonzero vector `b`.

# References: [Wikipedia](https://en.wikipedia.org/wiki/Vector_projection)

See also: [`sproj`](@ref), [`oproj`](@ref)
"""
proj(a, b) = (a ⋅ b / (b ⋅ b)) .* b

"""
Vector rejection
"""
oproj(a, b) = a .- (a ⋅ b / (b ⋅ b)) .* b

"""
Scalar projection
"""
sproj(a, b) = dot(a, b) / norm(b)
