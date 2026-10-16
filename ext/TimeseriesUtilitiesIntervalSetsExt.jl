module TimeseriesUtilitiesIntervalSetsExt

import TimeseriesUtilities as TU
using IntervalSets: AbstractInterval, leftendpoint, rightendpoint, isleftclosed, isrightclosed

function TU._range_indices(keys, i::AbstractInterval)
    t0, t1 = leftendpoint(i), rightendpoint(i)
    lo = isleftclosed(i) ? searchsortedfirst(keys, t0) : searchsortedlast(keys, t0) + 1
    hi = isrightclosed(i) ? searchsortedlast(keys, t1) : searchsortedfirst(keys, t1) - 1
    return lo:hi
end

TU.tmask!(A, i::AbstractInterval; dim = nothing) = TU.tmask!(A, (i,); dim)

end
