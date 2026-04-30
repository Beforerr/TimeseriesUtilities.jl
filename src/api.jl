# Common API for different types of arrays
# AxisKeys.jl: https://github.com/mcabbott/AxisKeys.jl
# DimensionalData.jl: https://github.com/rafaqz/DimensionalData.jl

"""
    dimnum(x, dim)

Get the ordinal of the dimension `dim` in `x`.
"""
dimnum(x, dim) = @something dim ndims(x)

function axiskeys end

function dims end

rebuild_axis(x, data, dim, keys) = data
rebuild_axis(x, dim, keys) = rebuild_axis(x, parent(x), dim, keys)

sorted_axis(sorted, dim; rev = false) = sorted

"""
    times(x)

Get time coordinate of `x`.
"""
times(x) = x

function samplingrate end
