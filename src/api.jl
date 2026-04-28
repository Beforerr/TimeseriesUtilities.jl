# Common API for different types of arrays
# AxisKeys.jl: https://github.com/mcabbott/AxisKeys.jl
# DimensionalData.jl: https://github.com/rafaqz/DimensionalData.jl

"""
    dimnum(x, dim)

Get the ordinal of the dimension `dim` in `x`.
"""
dimnum(x, dim) = @something dim ndims(x)

function set end

function axiskeys end

function dims end

rebuild_axis(x, data, dim, keys) = data

sorted_axis(sorted, dim; rev = false) = sorted

"""
    times(x)

Get time coordinate of `x`.
"""
times(x) = x

function samplingrate end

function _deriv_tfunc end

function tinterp_nans end

function workload_interp_setup end
