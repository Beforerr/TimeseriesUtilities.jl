# ============================================================================
# DimensionalData.jl Implementation of TimeseriesUtilities Interface
# ============================================================================

# --- Data Access ---

unwrap(x::AbstractDimArray) = parent(x)
unwrap(x::Dimension) = parent(lookup(x))

# --- Dimension Querying ---

dimnum(x::AbstractDimArray, query) = DimensionalData.dimnum(x, @something(query, TimeDim))
dimnum(x::AbstractDimStack, query) = DimensionalData.dimnum(x, @something(query, TimeDim))

dims(x::AbstractDimArray, dim) = DimensionalData.dims(x, dim)
dims(x::AbstractDimStack, dim) = DimensionalData.dims(x, dim)

axiskeys(x::AbstractDimArray, dim) = unwrap(DimensionalData.dims(x, dim))
axiskeys(x::AbstractDimStack, dim) = unwrap(DimensionalData.dims(x, dim))

# --- Array Reconstruction ---

rebuild(x::AbstractDimArray, data) = DimensionalData.rebuild(x; data = data)
rebuild(x::AbstractDimArray, data, newdims) = DimensionalData.rebuild(x, data, newdims)
rebuild(dim::Dimension, values) = DimensionalData.rebuild(dim, values)

set(x::AbstractDimArray, pair::Pair) = DimensionalData.set(x, pair)
set(x::AbstractDimStack, pair::Pair) = DimensionalData.set(x, pair)

"""
Returns the time indices of `x`.
"""
times(x::AbstractDimArray, args...) = lookup(timedim(x, args...))
times(x::AbstractDimStack, args...) = lookup(timedim(x, args...))

function timedim(x::Union{AbstractDimArray, AbstractDimStack}, query = nothing)
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

resolution(da::AbstractDimArray; kwargs...) = resolution(times(da); kwargs...)

# --- DimensionalData-specific Optimizations ---

# This is faster than `DimensionalData.format(rebuild(dim, x))` if the lookup trait keeps the same
fast_rebuild_dim(dim::Dimension, x) = rebuild(dim, rebuild(dim.val; data = x))

# --- DimensionalData-specific Overloads ---

function tstat(f, ds::AbstractDimStack, args...; query = nothing, dim = nothing)
    dim = @something dim dimnum(ds, query)
    return DimensionalData.maplayers(ds) do layer
        tstat(f, layer, args...; dim)
    end
end

function tinterp(A::AbstractDimArray, t; query = nothing, dim = nothing, kws...)
    dim = @something dim dimnum(A, query)
    out = tinterp(parent(A), unwrap(dims(A, dim)), t; dim, kws...)
    return if t isa AbstractTime
        out
    else
        newdim = fast_rebuild_dim(dims(A, dim), t)
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