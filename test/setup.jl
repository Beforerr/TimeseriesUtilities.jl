using Dates
using DimensionalData
using DimensionalData: Ti, Y

function workload_interp_setup(n = 4)
    times1 = DateTime(2020, 1, 1) + Day.(0:(n - 1))
    times2 = DateTime(2020, 1, 2) + Day.(0:(n - 1))
    times3 = DateTime(2020, 1, 1, 12) + Day.(0:(n - 2))

    da1 = DimArray(1:n, (Ti(times1),))
    da2 = DimArray(10:(10 + n - 1), (Ti(times2),))
    da3 = DimArray(hcat(5:(5 + n - 2), 8:2:(8 + 2n - 4)), (Ti(times3), Y([1, 2])))
    return da1, da2, da3
end
