module TimeseriesUtilitiesAxisKeysExt

import AxisKeys
import TimeseriesUtilities:
    axiskeys,
    dimnum,
    rebuild,
    rekey,
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

function rebuild(x::KeyedArray, data, dim::Integer, keys)
    newkeys = ntuple(ndims(x)) do i
        i == dim ? keys : axiskeys(x, i)
    end
    return _keyed(data, AxisKeys.dimnames(x), newkeys)
end

rebuild(x::KeyedArray, data) = _keyed(data, AxisKeys.dimnames(x), AxisKeys.axiskeys(x))
rekey(x::KeyedArray, dim, keys) = rebuild(x, unwrap(x), dim, keys)

rebuild(x::KeyedArray, data, dims, keys) = _keyed(data, map(d -> _name(x, d), dims), keys)

end
