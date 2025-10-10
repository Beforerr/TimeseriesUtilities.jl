"""
    tsubtract(x, f=nanmedian; dims=timedim(x))

Subtract a statistic (default function `f`: `nanmedian`) along dimensions (default: time dimension) from `x`.
"""
function tsubtract(x, f=nanmedian; dims=timedim(x))
    return x .- f(parent(x); dims=dimnum(x, dims))
end

function _deriv_tfunc(T1::Type, T2::Type)
    (T1 <: Real) && (T2 <: Dates.AbstractTime) ? Unitful.Quantity : identity
end

_deriv_tfunc(A, t) = _deriv_tfunc(eltype(A), eltype(t))

function _tderiv(A, times; dim=1)
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
    tderiv(A; dim = nothing, query = nothing)

Compute the time derivative of `A`. Set `lazy=true` for lazy evaluation.

See also: [deriv_data - PySPEDAS](https://pyspedas.readthedocs.io/en/latest/_modules/pyspedas/analysis/deriv_data.html)
"""
@inline function tderiv(args...; lazy=false, kw...)
    f = lazy ? DiffQ : _tderiv
    return f(args...; kw...)
end

function tderiv(A; dim = nothing, query = nothing, lazy=false)
    dim = @something dim dimnum(A, query)
    tdim = dims(A, dim)
    f = lazy ? DiffQ : _tderiv
    out = f(parent(A), unwrap(tdim); dim)
    newdims = ntuple(i -> i == dim ? @view(tdim[1:end-1]) : dims(A, i), ndims(A))
    return rebuild(A, out, newdims)
end


"""
    tnorm(A; dim=nothing, query=nothing)

Compute the norm of each slice in `A` along the specified dimension `dim` or `query`.

See also: [`tnorm_combine`](@ref)
"""
function tnorm(A; dim=nothing, query=nothing)
    dims = dimquery(dim, query)
    return norm.(eachslice(A; dims))
end

cross3(x, y) = cross(SV3(x), SV3(y))

"""
    tcross(x, y; dim = nothing, query=nothing)

Compute the cross product of two (arrays of) vectors along the specified dimension `dim` or `query`.

References:

  - https://docs.xarray.dev/en/stable/generated/xarray.cross.html
"""
function tcross(x, y; dim=nothing, query=nothing)
    dims = dimquery(dim, query)
    z = similar(x)
    map!(cross3, eachslice(z; dims), eachslice(x; dims), eachslice(y; dims))
    return z
end

"""
    tdot(x, y; dim=nothing, query=nothing)

Dot product of two arrays `x` and `y` along the `dim` dimension.
"""
function tdot(x, y; dim=nothing, query=nothing)
    dims = dimquery(dim, query)
    return dot.(eachslice(x; dims), eachslice(y; dims))
end

@inline function norm_combine(x::AbstractMatrix, dims)
    nn = norm.(eachslice(x; dims))
    return dims == 1 ?  hcat(x, nn) : vcat(x, nn')
end

"""
    tnorm_combine(x; dim=nothing, name=:magnitude)

Calculate the norm of each slice along `query` dimension and combine it with the original components.
"""
function tnorm_combine(x; dim=nothing, query=nothing, name=:magnitude)
    dim = @something dim dimnum(x, query)
    data = norm_combine(parent(x), dim)

    # Replace the original dimension with our new one that includes the magnitude
    odim = otherdims(x, dim) |> only
    odimType = basetypeof(odim)
    new_odim = odimType(vcat(odim.val, name))
    new_dims = map(d -> d isa odimType ? new_odim : d, dims(x))
    return rebuild(x, data, new_dims)
end


"""
    proj(a, b)

Vector projection of a vector `a` on (or onto) a nonzero vector `b`.

# References: [Wikipedia](https://en.wikipedia.org/wiki/Vector_projection)

See also: [`sproj`](@ref), [`oproj`](@ref)
"""
proj(a, b) = (a ⋅ b / (b ⋅ b)) .* b

function proj!(c, a, b)
    c .= (a ⋅ b / (b ⋅ b)) .* b
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
    tsproj(A, B; dim=nothing, query=nothing)

Compute scalar projection of `A` onto `B` along specified dimension `dim` or `query`.
"""
function tsproj(A, B; dim=nothing, query=nothing)
    dims = dimquery(dim, query)
    return sproj.(eachslice(A; dims), eachslice(B; dims))
end

"""
    tproj(A, B; dim=nothing, query=nothing)

Compute vector projection of `A` onto `B` along specified dimension `dim` or `query`.
"""
function tproj(A, B; dim=nothing, query=nothing)
    dims = dimnum(A, dimquery(dim, query))
    C = similar(A)
    as = eachslice(parent(A); dims)
    bs = eachslice(parent(B); dims)
    cs = eachslice(parent(C); dims)
    for i in eachindex(as, bs, cs)
        proj!(cs[i], as[i], bs[i])
    end
    return C
end


"""
    toproj(A, B; dim=nothing, query=nothing)

Compute vector rejection (orthogonal projection) of array `A` from `B` along specified dimension `dim` or `query`.
"""
function toproj(A, B; dim=nothing, query=nothing)
    dims = dimnum(A, dimquery(dim, query))
    C = similar(A)
    as = eachslice(parent(A); dims)
    bs = eachslice(parent(B); dims)
    cs = eachslice(parent(C); dims)
    for i in eachindex(as, bs, cs)
        oproj!(cs[i], as[i], bs[i])
    end
    return C
end

for f in (:tsubtract,)
    @eval $f(args...; kwargs...) = x -> $f(x, args...; kwargs...)
end
