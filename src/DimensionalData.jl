unwrap(x::AbstractDimArray) = parent(x)
unwrap(x::Dimension) = parent(lookup(x))

dimnum(x, query) = DimensionalData.dimnum(x, @something(query, TimeDim))

"""
Returns the time indices of `x`.
"""
times(x::AbstractDimArray, args...) = lookup(timedim(x, args...))
times(x::AbstractDimStack, args...) = lookup(timedim(x, args...))

function timedim(x, query = nothing)
    query = something(query, TimeDim)
    qdim = dims(x, query)
    return isnothing(qdim) ? dims(x, 1) : qdim
end

timerange(times::DimensionalData.Sampled) = timerange(parent(times))
timerange(times::Dimension) = timerange(parent(times))
timerange(x::AbstractDimArray) = timerange(times(x))

function groupby_dynamic(x::Dimension, args...; kwargs...)
    return groupby_dynamic(parent(lookup(x)), args...; kwargs...)
end

# This is faster than `DimensionalData.format(rebuild(dim, x))` if the lookup trait keeps the same
fast_rebuild_dim(dim, x) = rebuild(dim, rebuild(dim.val; data = x))

resolution(da::AbstractDimArray; kwargs...) = resolution(times(da); kwargs...)

function tinterp(A, t; query = nothing, dim = nothing, kws...)
    dim = @something dim dimnum(A, query)
    out = tinterp(parent(A), unwrap(dims(A, dim)), t; dim, kws...)
    return if t isa AbstractTime
        out
    else
        newdim = rebuild(dims(A, dim), t)
        newdims = ntuple(i -> i == dim ? newdim : dims(A, i), ndims(A))
        rebuild(A, out, DimensionalData.format(newdims, out))
    end
end

function tresample(A, dt; query = nothing, dim = nothing, kws...)
    dim = @something dim dimnum(A, query)
    old_times = unwrap(dims(A, dim))
    return tinterp(A, time_grid(old_times, dt); dim, kws...)
end

"""
    tinterp(A, B::AbstractDimArray; kws...)

Interpolate `A` to times in `B`
"""
tinterp(A, B::AbstractDimArray; kws...) = tinterp(A, times(B); kws...)