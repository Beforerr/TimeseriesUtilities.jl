unwrap(x::AbstractDimArray) = parent(x)
unwrap(x::Dimension) = parent(lookup(x))

dimnum(x, ::Nothing) = DimensionalData.dimnum(x, TimeDim)
dimnum(x, dim) = DimensionalData.dimnum(x, dim)

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

function smooth(da::AbstractDimArray, window; dim = nothing, kwargs...)
    dnum = dimnum(da, dim)
    tdim = DimensionalData.dims(da, dnum)
    data = smooth(parent(da), parent(tdim), window; dim = dnum, kwargs...)
    return rebuild(da; data)
end

function tinterp(A, t; dim = nothing, kws...)
    d = dimnum(A, dim)
    out = tinterp(parent(A), unwrap(dims(A, d)), t; dim = d, kws...)
    return if t isa AbstractTime
        out
    else
        newdim = rebuild(dims(A, d), t)
        newdims = ntuple(i -> i == d ? newdim : dims(A, i), ndims(A))
        rebuild(A, out, DimensionalData.format(newdims, out))
    end
end

function tresample(A, dt; dim = nothing, kws...)
    d = dimnum(A, dim)
    old_times = unwrap(dims(A, d))
    return tinterp(A, time_grid(old_times, dt); dim = d, kws...)
end

"""
    tinterp(A, B::AbstractDimArray; kws...)

Interpolate `A` to times in `B`
"""
tinterp(A, B::AbstractDimArray; kws...) = tinterp(A, times(B); kws...)

"""
    dropna(A::AbstractDimArray; dim=nothing)
    dropna(ds::DimStack; dim=nothing)

Remove slices containing NaN values along the `dim` dimension (defaults to the time dimension).
"""
function dropna(A::AbstractDimArray; dim = nothing)
    d = dimnum(A, dim)
    return _dropna(A; dim = d)
end

function dropna(ds::DimStack; dim = nothing)
    tdim = timedim(ds, dim)
    Dim, T = dimtype_eltype(tdim)
    odims = otherdims(ds, tdim)
    valid_idx = mapreduce(.*, values(ds)) do A
        vec(all(!isnan, A; dims = odims))
    end
    return ds[Dim(valid_idx)]
end

"""
    tinterp_nans(da::AbstractDimArray; dim=nothing, kwargs...)

Interpolate only the NaN values in `da` along the time dimension (or `dim`).
Non-NaN values are preserved exactly as they are.
"""
function tinterp_nans(da::AbstractDimArray; dim = nothing, kwargs...)
    u = parent(da)
    tdim = timedim(da, dim)
    t = parent(lookup(tdim))
    new_data = mapslices(u; dims = dimnum(da, dim)) do slice
        interpolate_nans(slice, t; kwargs...)
    end
    return rebuild(da; data = new_data)
end

for f in (:smooth, :tfilter)
    @eval $f(args...; kwargs...) = da -> $f(da, args...; kwargs...)
end

function tselect(da::AbstractDimArray, t; dim = nothing)
    Dim, T = dimtype_eltype(da, dim)
    return da[Dim(Near(T(t)))]
end
