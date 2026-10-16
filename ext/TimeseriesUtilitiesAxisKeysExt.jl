module TimeseriesUtilitiesAxisKeysExt

import AxisKeys
import TimeseriesUtilities:
    axiskeys,
    dimnum,
    rebuild,
    tdimnum,
    unwrap
using AxisKeys: KeyedArray
using TimeseriesUtilities: _axis, _axiskeys

unwrap(x::KeyedArray) = AxisKeys.keyless_unname(x)
dimnum(x::KeyedArray, name::Symbol) = AxisKeys.dim(x, name)
tdimnum(x::KeyedArray) = (i = AxisKeys.NamedDims.dim_noerror(AxisKeys.dimnames(x), :time); iszero(i) ? nothing : i)
axiskeys(x::KeyedArray, dim) = AxisKeys.axiskeys(x, dim)

_name(x, ::Val{i}) where {i} = AxisKeys.dimnames(x)[i]
_name(_, name::Symbol) = name
_name(_, ::Nothing) = :_

# `NamedDimsArray` accepts repeated `:_` names, keyword construction does not
@inline _keyed(data, names, keys) = KeyedArray(AxisKeys.NamedDims.NamedDimsArray(data, names), keys)

rebuild(x::KeyedArray, data, dims::Tuple) = _keyed(data, map(d -> _name(x, _axis(d)), dims), map(d -> _axiskeys(x, d), dims))

end
