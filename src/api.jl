# Common API for different types of arrays
# AxisKeys.jl: https://github.com/mcabbott/AxisKeys.jl
# DimensionalData.jl: https://github.com/rafaqz/DimensionalData.jl

import DimensionalData as DD

"""
    dimnum(x, query)

Get the number(s) of Dimension(s) as ordered in the dimensions of an object.

Extend the function for custom type `x`. By default, we fall back to `DimensionalData.dimnum`.
"""
function dimnum end

function set end

function axiskeys end

function dims end

function axiskeys(x::AbstractDimArray, dim)
    return unwrap(DD.dims(x, dim))
end

for f in (:set, :dims)
    @eval $f(args...; kwargs...) = DD.$f(args...; kwargs...)
end
