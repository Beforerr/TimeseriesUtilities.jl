"""
    tsplit(t0, t1, n::Int)
    tsplit(t0, t1, dt)

Split the range from `t0` to `t1` into `n` parts or `dt`-sized parts.
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
    n = ceil(Int, (t1 - t0) / dt)
    return tsplit(t0, t1, n)
end


"""
    unwrap(x)

Return the innermost object of the wrapped object `x` with similar behavior as `x` (e.g. same size, same type, etc.)
"""
unwrap(x) = x

dimquery(dim, query) = @something dim something(query, TimeDim)

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

using HybridArrays

rawview(x) = x

function rawview(x::AbstractArray{T}) where T
    return isbitstype(T) && sizeof(T) == sizeof(Int64) ? reinterpret(Int64, x) : x
end

rawview(x::AbstractTime) = Dates.value(x)

function hybridify(A, dims)
    sizes = ntuple(ndims(A)) do i
        i in dims ? StaticArrays.Dynamic() : size(A, i)
    end
    HybridArray{Tuple{sizes...}}(A)
end