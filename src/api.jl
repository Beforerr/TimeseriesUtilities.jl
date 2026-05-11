# Common API for different types of arrays
# AxisKeys.jl: https://github.com/mcabbott/AxisKeys.jl
# DimensionalData.jl: https://github.com/rafaqz/DimensionalData.jl

"""
    dimnum(x, dim)

Get the ordinal of the dimension `dim` in `x`.
"""
dimnum(x, dim) = @something dim ndims(x)

axiskeys(x, _) = throw(MethodError(axiskeys, (x,)))

dims(x, _) = throw(MethodError(dims, (x,)))

rebuild_axis(x, data, dim, keys) = data
rebuild_axis(x, dim, keys) = rebuild_axis(x, parent(x), dim, keys)
rebuild_axes(x, data, dims, keys) = data

sorted_axis(sorted, dim; rev = false) = sorted

"""
    times(x)

Get time coordinate of `x`.
"""
times(x) = x

function samplingrate end

_seconds(dt::Period) = dt / Second(1)
_seconds(dt) = dt

"""
    unwrap(x)

Return the innermost object (array) of wrapped `x` with similar behavior (e.g. same size, same type, etc.)
"""
unwrap(x) = x
