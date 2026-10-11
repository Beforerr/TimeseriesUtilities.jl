module TimeseriesUtilitiesAxisKeysExt

import AxisKeys
import TimeseriesUtilities:
    axiskeys,
    dimnum,
    rebuild,
    tdimnum,
    unwrap
using AxisKeys: KeyedArray

unwrap(x::KeyedArray) = AxisKeys.keyless_unname(x)
dimnum(x::KeyedArray, name::Symbol) = AxisKeys.dim(x, name)
tdimnum(x::KeyedArray) = (i = AxisKeys.NamedDims.dim_noerror(AxisKeys.dimnames(x), :time); iszero(i) ? nothing : i)
axiskeys(x::KeyedArray, dim) = AxisKeys.axiskeys(x, dim)

_name(x, ::Val{i}) where {i} = AxisKeys.dimnames(x)[i]
_name(_, name::Symbol) = name
_name(_, ::Nothing) = :_

# `NamedDimsArray` accepts repeated `:_` names, keyword construction does not
@inline _keyed(data, names, keys) = KeyedArray(AxisKeys.NamedDims.NamedDimsArray(data, names), keys)

@inline rebuild(x::KeyedArray, data, dim::Integer, keys) =
    _keyed(data, AxisKeys.dimnames(x), Base.setindex(AxisKeys.axiskeys(x), keys, dim))

rebuild(x::KeyedArray, data) = _keyed(data, AxisKeys.dimnames(x), AxisKeys.axiskeys(x))

rebuild(x::KeyedArray, data, dims, keys) = _keyed(data, map(d -> _name(x, d), dims), keys)

end
