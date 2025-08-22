# https://github.com/JuliaArrays/LazyArrays.jl/blob/master/src/lazyoperations.jl
struct Diff{T, N, Arr} <: AbstractArray{T, N}
    v::Arr
    dim::Int
end

function Diff(v)
    T = Base.eltype(v)
    DT = Base.promote_op(-, T, T)
    return Diff{DT, ndims(v), typeof(v)}(v, 1)
end

function Base.size(A::Diff)
    return size(A.v) .- ntuple(==(A.dim), ndims(A))
end

@propagate_inbounds function Base.getindex(A::Diff, i::Vararg{Int, N}) where {N}
    i_next = ntuple(==(A.dim), N) .+ i
    return A.v[i_next...] - A.v[i...]
end

# dx[i] = v[i+n] - v[i]
struct ForwardDiff{T, N, Arr} <: AbstractArray{T, N}
    v::Arr
    dim::Int
    n::Int
end

function ForwardDiff(v, dim, n)
    T = Base.eltype(v)
    DT = Base.promote_op(-, T, T)
    return ForwardDiff{DT, ndims(v), typeof(v)}(v, dim, n)
end

function Base.size(A::ForwardDiff)
    return size(A.v) .- A.n .* ntuple(==(A.dim), ndims(A))
end

@propagate_inbounds function Base.getindex(A::ForwardDiff, i::Vararg{Int, N}) where {N}
    i_next = A.n .* ntuple(==(A.dim), N) .+ i
    return A.v[i_next...] - A.v[i...]
end

# dx[i] = v[i+n] - v[i-n]
struct CentralDiff{T, N, Arr} <: AbstractArray{T, N}
    v::Arr
    dim::Int
    n::Int
end

function CentralDiff(v, dim, n)
    T = Base.eltype(v)
    DT = Base.promote_op(-, T, T)
    return CentralDiff{DT, ndims(v), typeof(v)}(v, dim, n)
end

function Base.size(A::CentralDiff)
    return size(A.v) .- 2 .* A.n .* ntuple(==(A.dim), ndims(A))
end

@propagate_inbounds function Base.getindex(A::CentralDiff, i::Vararg{Int, N}) where {N}
    i_next = 2 * A.n .* ntuple(==(A.dim), N) .+ i
    return A.v[i_next...] - A.v[i...]
end

"""
    DiffQ(v, t; dim=1)

Difference quotient of `v` with respect to `t`.

To avoid undefined behavior for division by Date/DateTime, we convert the time difference to a `Unitful.Quantity` if `eltype(v)` is not a `Unitful.Quantity`.
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
