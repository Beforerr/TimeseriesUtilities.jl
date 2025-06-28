# https://github.com/JuliaArrays/LazyArrays.jl/blob/master/src/lazyoperations.jl
struct Diff{T,N,Arr} <: AbstractArray{T,N}
    v::Arr
    dim::Int
end

function Diff(v)
    T = Base.eltype(v)
    DT = Base.promote_op(-, T, T)
    return Diff{DT,ndims(v),typeof(v)}(v, 1)
end

function Base.size(A::Diff)
    return size(A.v) .- ntuple(==(A.dim), ndims(A))
end

function Base.getindex(A::Diff, i::Vararg{Int,N}) where {N}
    i_next = ntuple(==(A.dim), N) .+ i
    return A.v[i_next...] - A.v[i...]
end

"""
    DiffQ(v, t; dims=1)

Difference quotient of `v` with respect to `t`.
"""
struct DiffQ{T,N,D,A1,A2} <: AbstractArray{T,N}
    v::A1
    t::A2
end

function DiffQ(v::AbstractArray{T1,N}, t::AbstractVector{T2}; dims=1) where {T1,N,T2}
    size(v, dims) == length(t) || throw(ArgumentError("v and t must have the same size"))
    _T1 = Base.promote_op(-, T1, T1)
    _T2 = Base.promote_op(-, T2, T2)
    T = Base.promote_op(/, _T1, _T2)
    return DiffQ{T,N,dims,typeof(v),typeof(t)}(v, t)
end

function Base.getindex(A::DiffQ{<:Any,N,D}, i::Vararg{Int,N}) where {N,D}
    i_next = ntuple(==(D), N) .+ i
    return (A.v[i_next...] - A.v[i...]) / (A.t[i_next[D]] - A.t[i[D]])
end

function Base.size(A::DiffQ{<:Any,N,D}) where {N,D}
    return size(A.v) .- ntuple(==(D), ndims(A))
end

tderiv(::Type{DiffQ}, data, times; dims=1) = DiffQ(data, times; dims)