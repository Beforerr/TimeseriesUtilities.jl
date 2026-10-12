"""
    tsubtract(f, x, [window]; dim=nothing)

Subtract from each channel of `x` its statistic `f` along time (or `dim`) over all samples, or with `window`.
"""
function tsubtract(f, x, window = nothing; dim = nothing)
    data = unwrap(x)
    isnothing(window) || return rebuild(x, data .- unwrap(smooth(f, x, window; dim)))
    d = resolve_dim(x, dim)
    stats = _map(f, eachslice(data; dims = other_dims(data, d)))
    return rebuild(x, data .- reshape(stats, Base.setindex(size(data), 1, d)))
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
    tfunc = _deriv_tfunc(A, times)
    # times vary along `dim` only
    shape = ntuple(i -> i == dim ? length(r0[dim]) : 1, N)
    t0, t1 = reshape(view(times, r0[dim]), shape), reshape(view(times, r1[dim]), shape)
    return (view(A, r1...) .- view(A, r0...)) ./ tfunc.(t1 .- t0)
end

"""
    tderiv(A; dim=nothing, lazy=false)

Forward-difference time derivative of `A`, keyed by the left time of each pair;
`Date`/`DateTime` differences are taken in seconds. `lazy=true` returns a lazy array.
"""
@inline function tderiv(A; lazy = false, dim = nothing)
    # not `_ontime`: `lazy` would not fold inside closure, so result type would be a `Union`
    d = resolve_dim(A, dim)
    f = lazy ? DiffQ : _tderiv
    t = axiskeys(A, d)
    return rebuild(A, f(unwrap(A), t; dim = d), d, @view t[1:(end - 1)])
end

# Vector projection, rejection and scalar projection of `a` on `b`.
proj(a, b) = (a ⋅ b / (b ⋅ b)) .* b
oproj(a, b) = a .- (a ⋅ b / (b ⋅ b)) .* b
sproj(a, b) = dot(a, b) / norm(b)

# https://github.com/JuliaArrays/LazyArrays.jl/blob/master/src/lazyoperations.jl
# https://gist.github.com/Beforerr/47218b68c975bdf68807ccbfe969f641
"""
    DiffQ(v, t; dim=1)

Difference quotient of `v` with respect to `t`.

Date/DateTime differences are converted to seconds before division.
"""
struct DiffQ{T, N, D, A1, A2} <: AbstractArray{T, N}
    v::A1
    t::A2
end

function DiffQ(v::AbstractArray{T1, N}, t::AbstractVector{T2}; dim = 1) where {T1, N, T2}
    size(v, dim) == length(t) || throw(ArgumentError("v and t must have the same size"))
    _T1 = Base.promote_op(-, T1, T1)
    _T2 = Base.promote_op(-, T2, T2)
    tfunc = _deriv_tfunc(v, t)
    T = Base.promote_op(/, _T1, Base.promote_op(tfunc, _T2))
    return DiffQ{T, N, dim, typeof(v), typeof(t)}(v, t)
end

@propagate_inbounds function Base.getindex(A::DiffQ{T, N, D}, i::Vararg{Int, N}) where {T, N, D}
    i_next = ntuple(==(D), N) .+ i
    tfunc = _deriv_tfunc(A.v, A.t)
    return (A.v[i_next...] - A.v[i...]) / tfunc(A.t[i_next[D]] - A.t[i[D]])
end

function Base.size(A::DiffQ{<:Any, N, D}) where {N, D}
    return size(A.v) .- ntuple(==(D), ndims(A))
end
