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
    dims,
    lookup,
    maplayers,
    otherdims,
    rebuild
import TimeseriesUtilities as TU
import TimeseriesUtilities:
    axiskeys,
    dimnum,
    dropna,
    groupby_dynamic,
    norm_combine,
    smooth,
    sorted_axis,
    tnorm_combine,
    tdimnum,
    tstat,
    unwrap

const AbstractDimLike = Union{AbstractDimArray, AbstractDimStack}

unwrap(x::AbstractDimArray) = parent(x)
unwrap(x::Dimension) = parent(lookup(x))
dimnum(x::AbstractDimLike, dim::DD.Dimensions.DimOrDimType) = DD.dimnum(x, dim)
_dimnum(x, dim) = DD.hasdim(x, dim) ? DD.dimnum(x, dim) : nothing
tdimnum(x::AbstractDimLike) = @something _dimnum(x, TimeDim) _dimnum(x, DD.Dim{:time}) Some(nothing)
@inline axiskeys(x::AbstractDimLike, dim) = unwrap(DD.dims(x, dim))

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
_newdim(x, ::Val{i}, keys) where {i} = _newdim(x, dims(x)[i], keys)
_newdim(_, dim::Dimension, keys) = rebuild(dim, _rekey(lookup(dim), keys))
_newdim(_, name::Symbol, keys) = DD.format(DD.Dim{name}(keys), axes(keys, 1))
_newdim(_, ::Nothing, keys) = DD.format(DD.AnonDim(keys), axes(keys, 1))

@inline TU.rebuild(x::AbstractDimArray, data, dim::Integer, keys) =
    rebuild(x, data, Base.setindex(dims(x), _newdim(x, dims(x)[dim], keys), dim))
TU.rebuild(x::AbstractDimArray, data, newdims, keys) =
    rebuild(x, data, map((d, k) -> _newdim(x, d, k), newdims, keys))

TU.rebuild(x::AbstractDimArray, data) = rebuild(x, data)
TU.rekey(x::AbstractDimArray, dim, keys) = TU.rebuild(x, parent(x), dim, keys)

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
