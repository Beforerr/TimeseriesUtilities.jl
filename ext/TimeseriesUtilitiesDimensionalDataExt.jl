module TimeseriesUtilitiesDimensionalDataExt

import DimensionalData as DD
using Dates
using DimensionalData:
    AbstractDimArray,
    AbstractDimStack,
    DimArray,
    Dimension,
    TimeDim,
    basetypeof,
    lookup,
    maplayers,
    otherdims,
    rebuild
import TimeseriesUtilities as TU
import TimeseriesUtilities:
    axiskeys,
    dimnum,
    dims,
    dropna,
    groupby_dynamic,
    norm_combine,
    smooth,
    sorted_axis,
    times,
    tnorm_combine,
    tstat,
    unwrap

const AbstractDimLike = Union{AbstractDimArray, AbstractDimStack}

unwrap(x::AbstractDimArray) = parent(x)
unwrap(x::Dimension) = parent(lookup(x))
dimnum(x::AbstractDimLike, dim) = DD.dimnum(x, @something dim TimeDim)
axiskeys(x::AbstractDimLike, dim) = unwrap(DD.dims(x, dim))
dims(x::AbstractDimLike, dim) = DD.dims(x, dim)
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

_span(keys::AbstractRange) = DD.Regular(step(keys))
_span(keys) = DD.Irregular((nothing, nothing))
_rekey(l::DD.Sampled, keys) = rebuild(l; data = keys, span = _span(keys))
_rekey(l, keys) = rebuild(l; data = keys)
_newdim(dim::Dimension, keys) = rebuild(dim, _rekey(lookup(dim), keys))
_newdim(dim::Symbol, keys) = DD.format(DD.Dim{dim}(keys), axes(keys, 1))

@inline TU.rebuild(x::AbstractDimArray, data, dim::Integer, keys) =
    rebuild(x, data, Base.setindex(DD.dims(x), _newdim(DD.dims(x, dim), keys), dim))
TU.rebuild(x::AbstractDimArray, data, newdims, keys) = rebuild(x, data, map(_newdim, newdims, keys))

TU.rebuild(x::AbstractDimArray, data) = rebuild(x, data)

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

function tstat(f, ds::AbstractDimStack, args...; dim = nothing)
    d = dimnum(ds, dim)
    return maplayers(ds) do layer
        tstat(f, layer, args...; dim = d)
    end
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

end
