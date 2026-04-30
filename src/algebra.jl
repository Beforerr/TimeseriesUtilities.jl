"""
    tsubtract(x, op=nanmedian; dim=nothing)

Subtract a statistic `op` along time dimension (or `dim`) from `x`.
"""
function tsubtract(x, op = nanmedian; dim = nothing)
    dims = dimnum(x, dim)
    return x .- op(parent(x); dims)
end

_seconds(dt::Period) = dt / Second(1)
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
    out = f(parent(A), times; dim)
    return rebuild_axis(A, out, dim, @view times[1:end-1])
end

"""
    tnorm(A; dim=nothing)

Compute the norm of each slice in `A` along dimension `dim`.

See also: [`tnorm_combine`](@ref)
"""
tnorm(A; dim = nothing) =
    norm.(eachslice(A; dims = dimnum(A, dim)))

cross3(x, y) = cross(SV3(x), SV3(y))

"""
    tcross(x, y; dim=nothing)

Compute the cross product of two (arrays of) vectors along dimension `dim`.
"""
function tcross(x, y; dim = nothing)
    dims = dimnum(x, dim)
    z = similar(x)
    map!(cross3, eachslice(z; dims), eachslice(x; dims), eachslice(y; dims))
    return z
end

"""
    tdot(x, y; dim=nothing)

Dot product of two arrays `x` and `y` along dimension `dim`.
"""
function tdot(x, y; dim = nothing)
    dims = dimnum(x, dim)
    return dot.(eachslice(x; dims), eachslice(y; dims))
end

@inline function norm_combine(x::AbstractMatrix, dims)
    nn = norm.(eachslice(x; dims))
    return dims == 1 ? hcat(x, nn) : vcat(x, nn')
end

"""
    tnorm_combine(x; dim=nothing, name=:magnitude)

Calculate the norm of each slice along `dim` and combine it with the original components.
"""
function tnorm_combine(x::AbstractMatrix; dim = nothing, name = :magnitude)
    d = dimnum(x, dim)
    return norm_combine(x, d)
end


"""
    proj(a, b)

Vector projection of a vector `a` on (or onto) a nonzero vector `b`.

# References: [Wikipedia](https://en.wikipedia.org/wiki/Vector_projection)

See also: [`sproj`](@ref), [`oproj`](@ref)
"""
proj(a, b) = (a ⋅ b / (b ⋅ b)) .* b

function proj!(c, a, b)
    return c .= (a ⋅ b / (b ⋅ b)) .* b
end

"""
Vector rejection
"""
oproj(a, b) = a .- (a ⋅ b / (b ⋅ b)) .* b

function oproj!(c, a, b)
    return c .= a .- (a ⋅ b / (b ⋅ b)) .* b
end

"""
Scalar projection
"""
sproj(a, b) = dot(a, b) / norm(b)

"""
    tsproj(A, B; dim=nothing)

Compute scalar projection of `A` onto `B` along dimension `dim`.
"""
function tsproj(A, B; dim = nothing)
    dims = dimnum(A, dim)
    return sproj.(eachslice(A; dims), eachslice(B; dims))
end

function _tproj(f!, A, B; dim = nothing)
    C = similar(A)
    dims = dimnum(A, dim)
    as, bs, cs = eachslice(A; dims), eachslice(B; dims), eachslice(C; dims)
    for i in eachindex(as, bs, cs)
        f!(cs[i], as[i], bs[i])
    end
    return C
end

"""
    tproj(A, B; dim=nothing)

Compute vector projection of `A` onto `B` along dimension `dim`.
"""
tproj(A, B; dim = nothing) = _tproj(proj!, A, B; dim)

"""
    toproj(A, B; dim=nothing)

Compute vector rejection (orthogonal projection) of `A` from `B` along dimension `dim`.
"""
toproj(A, B; dim = nothing) = _tproj(oproj!, A, B; dim)
