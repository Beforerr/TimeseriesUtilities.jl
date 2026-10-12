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
    otherdims
import TimeseriesUtilities:
    axiskeys,
    dimnum,
    dropna,
    rebuild,
    resolve_dim,
    sorted_axis,
    tdimnum,
    tstat,
    unwrap

const AbstractDimLike = Union{AbstractDimArray, AbstractDimStack}

unwrap(x::AbstractDimArray) = parent(x) 
unwrap(x::Dimension) = (l = lookup(x); l isa DD.Lookup ? parent(l) : l)
dimnum(x::AbstractDimLike, dim::DD.Dimensions.DimOrDimType) = DD.dimnum(x, dim)
_dimnum(x, dim) = DD.hasdim(x, dim) ? DD.dimnum(x, dim) : nothing
tdimnum(x::AbstractDimLike) = @something _dimnum(x, TimeDim) _dimnum(x, DD.Dim{:time}) Some(nothing)
@inline axiskeys(x::AbstractDimLike, dim) = parent(lookup(DD.dims(x, dim)))
# `DD.dims(x, ::Int)` does not fold a constant index when inferred nested
@inline axiskeys(x::AbstractDimLike, dim::Integer) = parent(lookup(dims(x)[dim]))

_ordered_lookup(lookup::DD.Sampled, ::Val{false}) = DD.rebuild(lookup; order = DD.ForwardOrdered())
_ordered_lookup(lookup::DD.Sampled, ::Val{true}) = DD.rebuild(lookup; order = DD.ReverseOrdered())
_ordered_lookup(lookup, rev) = lookup

@inline function sorted_axis(sorted::AbstractDimArray, D; rev = false)
    olddims = DD.dims(sorted)
    olddim = olddims[D]
    newdim = DD.rebuild(olddim, _ordered_lookup(lookup(olddim), Val(rev)))
    newdims = Base.setindex(olddims, newdim, D)
    return DD.rebuild(sorted, parent(sorted), newdims)
end

@inline function sorted_axis(sorted::AbstractDimStack, D; rev = false)
    olddims = DD.dims(sorted)
    olddim = olddims[D]
    newdim = DD.rebuild(olddim; val = _ordered_lookup(lookup(olddim), Val(rev)))
    newdims = Base.setindex(olddims, newdim, D)
    return DD.rebuild(sorted; dims = newdims)
end

_span(keys::AbstractRange) = DD.Regular(step(keys))
_span(keys) = DD.Irregular((nothing, nothing))
_rekey(l::DD.Sampled, keys) = DD.rebuild(l; data = keys, span = _span(keys))
_rekey(l, keys) = DD.rebuild(l; data = keys)
_newdim(x, ::Val{i}, keys) where {i} = _newdim(x, dims(x)[i], keys)
_newdim(_, dim::Dimension, keys) = DD.rebuild(dim, _rekey(lookup(dim), keys))
_newdim(_, name::Symbol, keys) = DD.format(DD.Dim{name}(keys), axes(keys, 1))
_newdim(_, ::Nothing, keys) = DD.format(DD.AnonDim(keys), axes(keys, 1))

@inline rebuild(x::AbstractDimArray, data, dim::Integer, keys) =
    DD.rebuild(x, data, Base.setindex(dims(x), _newdim(x, dims(x)[dim], keys), dim))
rebuild(x::AbstractDimArray, data, newdims, keys) =
    DD.rebuild(x, data, map((d, k) -> _newdim(x, d, k), newdims, keys))

rebuild(x::AbstractDimArray, data) = DD.rebuild(x, data)

"""
    dropna(ds::AbstractDimStack; dim=nothing)

Remove slices containing NaN values along the `dim` dimension.
"""
function dropna(ds::AbstractDimStack; dim = nothing)
    d = resolve_dim(ds, dim)
    tdim = DD.dims(ds, d)
    odims = otherdims(ds, d)
    valid_idx = mapreduce(.*, values(ds)) do A
        vec(all(!isnan, A; dims = odims))
    end
    return ds[basetypeof(tdim)(valid_idx)]
end

function tstat(f, ds::AbstractDimStack, args...; dim = nothing, kw...)
    d = resolve_dim(ds, dim)
    return maplayers(ds) do layer
        tstat(f, layer, args...; dim = d, kw...)
    end
end

end
