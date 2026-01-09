module TimeseriesUtilitiesAxisKeysExt

import TimeseriesUtilities as TU
import AxisKeys
import TimeseriesUtilities: unwrap, dimnum, set, dims, axiskeys
using AxisKeys: KeyedArray

TU.unwrap(x::KeyedArray) = AxisKeys.keyless_unname(x)
TU.dimnum(x::KeyedArray, query) = AxisKeys.dim(x, @something(query, :time))
TU.axiskeys(x::KeyedArray, dim) = AxisKeys.axiskeys(x, dim)
TU.dims(x::KeyedArray, dim) = AxisKeys.dimnames(x, dim)

function TU.set(x::KeyedArray, pair::Pair)
    dim, new_keys = pair
    dn = dimnum(x, dim)
    new_axiskeys = ntuple(ndims(x)) do i
        i == dn ? new_keys : axiskeys(x, i)
    end
    return KeyedArray(parent(x), new_axiskeys)
end

end
