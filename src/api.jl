# Common API for different types of arrays
# AxisKeys.jl: https://github.com/mcabbott/AxisKeys.jl
# DimensionalData.jl: https://github.com/rafaqz/DimensionalData.jl

"""
    dimnum(x, query)

Get the number(s) of Dimension(s) as ordered in the dimensions of an object.

Extend the function for custom type `x`.
"""
function dimnum end

function set end

function axiskeys end

function dims end

"""
    times(x)

Get the time indices of `x`. Extend for custom types.
"""
function times end

"""
    timedim(x, query=nothing)

Get the time dimension of `x`. Extend for custom types.
"""
function timedim end

"""
    rebuild(x, data)
    rebuild(x, data, dims)

Rebuild `x` with new `data` (and optionally new `dims`).
Generic fallback returns `data` as-is. Override for custom types to preserve metadata.
"""
rebuild(x, data) = data
rebuild(x, data, dims) = data
