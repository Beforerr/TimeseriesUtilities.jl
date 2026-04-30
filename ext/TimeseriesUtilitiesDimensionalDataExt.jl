module TimeseriesUtilitiesDimensionalDataExt

import DimensionalData as DD
using Dates
using DimensionalData:
    AbstractDimArray,
    AbstractDimStack,
    DimArray,
    Dimension,
    Ti,
    TimeDim,
    Y,
    basetypeof,
    lookup,
    maplayers,
    otherdims,
    rebuild
import TimeseriesUtilities as TU
import TimeseriesUtilities:
    DiffQ,
    _tderiv,
    axiskeys,
    dimnum,
    dims,
    dropna,
    find_outliers,
    groupby_dynamic,
    norm_combine,
    rebuild_axis,
    resolution,
    smooth,
    sorted_axis,
    tderiv,
    tfilter,
    times,
    tnorm_combine,
    tstat,
    unwrap,
    workload_interp_setup

const AbstractDimLike = Union{AbstractDimArray, AbstractDimStack}

unwrap(x::AbstractDimArray) = parent(x)
unwrap(x::Dimension) = parent(lookup(x))
dimnum(x::AbstractDimLike, dim) = DD.dimnum(x, @something dim TimeDim)
axiskeys(x::AbstractDimLike, dim) = unwrap(DD.dims(x, dim))
times(x::AbstractDimLike, dim = nothing) = axiskeys(x, dimnum(x, dim))

for f in (:dims,)
    @eval $f(args...; kwargs...) = DD.$f(args...; kwargs...)
end

function groupby_dynamic(x::Dimension, args...; kwargs...)
    return groupby_dynamic(parent(lookup(x)), args...; kwargs...)
end

_ordered_lookup(lookup::DD.Sampled, ::Val{false}) = rebuild(lookup; order = DD.ForwardOrdered())
_ordered_lookup(lookup::DD.Sampled, ::Val{true}) = rebuild(lookup; order = DD.ReverseOrdered())
_ordered_lookup(lookup, rev) = lookup

@inline function sorted_axis(sorted::AbstractDimArray, D; rev = false)
    olddims = DD.dims(sorted)
    olddim = olddims[D]
    newdim = rebuild(olddim, _ordered_lookup(lookup(olddim), Val(rev)))
    newdims = Base.setindex(olddims, newdim, D)
    return rebuild(sorted, parent(sorted), newdims)
end

@inline function sorted_axis(sorted::AbstractDimStack, D; rev = false)
    olddims = DD.dims(sorted)
    olddim = olddims[D]
    newdim = rebuild(olddim; val = _ordered_lookup(lookup(olddim), Val(rev)))
    newdims = Base.setindex(olddims, newdim, D)
    return rebuild(sorted; dims = newdims)
end

@inline function rebuild_axis(x::AbstractDimArray, data, dim, keys)
    olddims = DD.dims(x)
    olddim = olddims[dim]
    newdim = DD.rebuild(olddim; val = DD.rebuild(lookup(olddim); data = keys))
    newdims = Base.setindex(olddims, newdim, dim)
    return rebuild(x, data, newdims)
end

resolution(da::AbstractDimArray; kwargs...) = resolution(times(da); kwargs...)

function smooth(da::AbstractDimArray, window; dim = nothing, kwargs...)
    dnum = dimnum(da, dim)
    tdim = DD.dims(da, dnum)
    data = smooth(parent(da), parent(tdim), window; dim = dnum, kwargs...)
    return rebuild(da; data)
end

"""
    dropna(ds::AbstractDimStack; dim=nothing)

Remove slices containing NaN values along the `dim` dimension.
"""
function dropna(ds::AbstractDimStack; dim = nothing)
    d = dimnum(ds, dim)
    tdim = DD.dims(ds, d)
    odims = otherdims(ds, d)
    valid_idx = mapreduce(.*, values(ds)) do A
        vec(all(!isnan, A; dims = odims))
    end
    return ds[basetypeof(tdim)(valid_idx)]
end

for f in (:smooth, :tfilter)
    @eval $f(args...; kwargs...) = da -> $f(da, args...; kwargs...)
end

function tstat(f, ds::AbstractDimStack, args...; dim = nothing)
    d = dimnum(ds, dim)
    return maplayers(ds) do layer
        tstat(f, layer, args...; dim = d)
    end
end

function tderiv(A::AbstractDimArray; dim = nothing, lazy = false)
    d = dimnum(A, dim)
    tdim = dims(A, d)
    f = lazy ? DiffQ : _tderiv
    out = f(parent(A), unwrap(tdim); dim = d)
    newdims = ntuple(i -> i == d ? @view(tdim[1:(end - 1)]) : dims(A, i), ndims(A))
    return rebuild(A, out, newdims)
end

function tnorm_combine(x::AbstractDimArray; dim = nothing, name = :magnitude)
    d = dimnum(x, dim)
    data = norm_combine(parent(x), d)

    odim = otherdims(x, d) |> only
    odimType = basetypeof(odim)
    new_odim = odimType(vcat(odim.val, name))
    new_dims = map(dd -> dd isa odimType ? new_odim : dd, dims(x))
    return rebuild(x, data, new_dims)
end

function find_outliers(A::AbstractDimArray, args...; dim = nothing, kw...)
    d = dimnum(A, dim)
    return find_outliers(parent(A), args...; dim = d, kw...)
end

function workload_interp_setup(n = 4)
    times1 = DateTime(2020, 1, 1) + Day.(0:(n - 1))
    times2 = DateTime(2020, 1, 2) + Day.(0:(n - 1))
    times3 = DateTime(2020, 1, 1, 12) + Day.(0:(n - 2))

    da1 = DimArray(1:n, (Ti(times1),))
    da2 = DimArray(10:(10 + n - 1), (Ti(times2),))
    da3 = DimArray(hcat(5:(5 + n - 2), 8:2:(8 + 2n - 4)), (Ti(times3), Y([1, 2])))
    return da1, da2, da3
end

end
