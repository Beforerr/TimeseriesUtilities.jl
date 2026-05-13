/ₜ(x, n) = x / n
/ₜ(x::P, n) where {P <: Dates.AbstractTime} = P(cld(Dates.value(x), n))

# Handle numeric offsets for datetime-like types (default to Unix epoch)
struct TimeOffsets{T, A <: AbstractArray} <: AbstractVector{T}
    offsets::A
    t0::T
end

TimeOffsets(offsets) = TimeOffsets(offsets, Dates.unix2datetime(0))

Base.size(to::TimeOffsets) = size(to.offsets)
function Base.getindex(to::TimeOffsets, i::Int)
    _add(t0, dt) = t0 + dt
    _add(t0::Dates.AbstractTime, dt::Number) = t0 + Nanosecond(round(Int, 1.0e9 * dt))
    return _add(to.t0, to.offsets[i])
end

"""
    tsplit(t0, t1, n::Int)
    tsplit(t0, t1, dt)
    tsplit(t0, t1, dtType::Type{<:Period})

Split the range from `t0` to `t1` into `n` parts, `dt`-sized parts, or by period type (e.g., Month, Day).
"""
tsplit(t0, t1, dt) = collect(IntervalRange(t0, t1, dt))
tsplit((t0, t1), arg) = tsplit(t0, t1, arg)

function stat_relerr(itr, f)
    m = f(itr)
    relerrs = abs.(extrema(itr) .- m) ./ m
    relerr = maximum(relerrs)
    return m, relerr
end

"""
    window_bf_sizes(window)

Converts a window specification to backward and forward window sizes.

When window is a positive integer scalar, the window is centered about the current element and contains window-1 neighboring elements.
If window is even, then the window is centered about the current and previous elements.
"""
function window_bf_sizes(window::Integer)
    return isodd(window) ? (window ÷ 2, window ÷ 2) : (window ÷ 2, window ÷ 2 - 1)
end

function window_bf_sizes(window)
    @assert length(window) == 2 "Window must be of length 2"
    return window
end

other_dims(A, dim) = filter(!=(dim), ntuple(identity, ndims(A)))

# Type-stable alternative to selectdim: ntuple with Val(N) makes every index
# position a compile-time constant, so the SubArray type is fully inferred.
@inline _vdim(A::AbstractArray{T, N}, ::Val{d}, k) where {T, N, d} =
    @inbounds view(A, ntuple(j -> j == d ? k : Colon(), Val(N))...)

# Lazy slice iterator for interpolators
# T encodes the slice representation — copy-based by default, SArray when
# the StaticArrays extension is loaded. Splines require T to support round-trip
# arithmetic (e.g. u[i+1]-u[i] returns the same type), so views are not suitable.
struct LazySlices{T, A, d} <: AbstractVector{T}
    data::A
end

function LazySlices(A::AbstractArray{Tv, N}, d) where {Tv, N}
    return LazySlices{Array{Tv, N - 1}, typeof(A), d}(A)
end

Base.length(s::LazySlices{T, A, d}) where {T, A, d} = size(s.data, d)

@inline function Base.getindex(s::LazySlices{T, A, d}, k::Int) where {T <: Array, A, d}
    return copy(_vdim(s.data, Val(d), k))
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
