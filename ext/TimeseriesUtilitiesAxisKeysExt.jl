module TimeseriesUtilitiesAxisKeysExt

import AxisKeys
import TimeseriesUtilities:
    axiskeys,
    dimnum,
    dims,
    rebuild_axis,
    times,
    unwrap
using AxisKeys: KeyedArray

unwrap(x::KeyedArray) = AxisKeys.keyless_unname(x)
dimnum(x::KeyedArray, dim) = AxisKeys.dim(x, @something dim :time)
axiskeys(x::KeyedArray, dim) = AxisKeys.axiskeys(x, dim)
dims(x::KeyedArray, dim) = AxisKeys.dimnames(x, dim)
times(x::KeyedArray, dim = nothing) = axiskeys(x, dimnum(x, dim))
function rebuild_axis(x::KeyedArray, data, dim, keys)
    names = ntuple(i -> dims(x, i), ndims(x))
    newkeys = ntuple(ndims(x)) do i
        i == dim ? keys : axiskeys(x, i)
    end
    return KeyedArray(data; NamedTuple{names}(newkeys)...)
end

end
