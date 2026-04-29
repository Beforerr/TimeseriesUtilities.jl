"""
    tsplit(t0, t1, n::Int)
    tsplit(t0, t1, dt)
    tsplit(t0, t1, dtType::Type{<:Period})

Split the range from `t0` to `t1` into `n` parts, `dt`-sized parts, or by period type (e.g., Month, Day).
"""
function tsplit(t0, t1, n::Int)
    return if n <= 1
        [(t0, t1)]
    else
        dt = (t1 - t0) / n
        map(1:n) do i
            t0 + (i - 1) * dt, min(t0 + i * dt, t1)
        end
    end
end

function tsplit(t0, t1, dt)
    T = promote_type(typeof(t0), typeof(t1))
    periods = NTuple{2, T}[]
    current = t0

    while current < t1
        next_period = current + dt
        push!(periods, (current, min(next_period, t1)))
        current = next_period
    end

    return periods
end

tsplit(t0, t1, dtType::Type{<:Period}) = tsplit(t0, t1, dtType(1))
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

# Lazy iterator that yields SArray slices along `dim` without heap-allocating views.
# ntuple(Val(n)) unrolls at compile time, so all intermediates are stack-allocated.
# Inspired by HybridArrays.jl
struct LazyStaticSlices{T, N, S <: SArray, dim, A <: AbstractArray{T, N}} <: AbstractVector{S}
    data::A
end

function LazyStaticSlices(A::AbstractArray{T, N}, dim::Int) where {T, N}
    other = ntuple(i -> size(A, i < dim ? i : i + 1), N - 1)
    S = SArray{Tuple{other...}, T, N - 1, prod(other)}
    return LazyStaticSlices{T, N, S, dim, typeof(A)}(A)
end

Base.length(s::LazyStaticSlices{T, N, S, d}) where {T, N, S, d} = size(s.data, d)

@inline function Base.getindex(s::LazyStaticSlices{T, N, S, d}, k::Int) where {T, N, S, d}
    slice_sz = size(S)
    return S(ntuple(Val(length(S))) do p
        ci = CartesianIndices(slice_sz)[p]
        full_idx = ntuple(Val(N)) do j
            j < d ? ci[j] : j == d ? k : ci[j - 1]
        end
        @inbounds s.data[full_idx...]
    end)
end
