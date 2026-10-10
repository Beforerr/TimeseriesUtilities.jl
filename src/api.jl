# Common API for different types of arrays
# AxisKeys.jl: https://github.com/mcabbott/AxisKeys.jl
# DimensionalData.jl: https://github.com/rafaqz/DimensionalData.jl

"""
    tdimnum(x)

Index of the time dimension of `x`, or `nothing` if `x` has none.
Arrays without a named time dimension take time along their last dimension.
"""
tdimnum(x) = nothing
tdimnum(x::AbstractArray) = ndims(x)

"""
    dimnum(x, dim)

Index of dimension `dim` in `x`; `nothing` selects the time dimension ([`tdimnum`](@ref)).

Containers with named dimensions add methods for their name types.
"""
dimnum(_, dim::Integer) = dim
dimnum(x, ::Nothing) =
    @something tdimnum(x) throw(ArgumentError("$(typeof(x)) has no time dimension; pass `dim`"))

function axiskeys end

# Containers without a `rebuild` method degrade to their data, dropping their axis keys.
rebuild(x, data) = data
# rebuild multiple dimensions; each entry of `dims` is a `Val{i}` (dimension `i` of `x`, static so
# the result type is inferable), a `Symbol` (new named dimension) or `nothing` (new anonymous dimension)
rebuild(x, data, dims, keys) = data
# rebuild one dimension with new data and keys
rebuild(x, data, dim::Integer, key) = data

# Replace the keys of dimension `dim`
rekey(x, dim, keys) = throw(ArgumentError("$(typeof(x)) has no axis keys to replace"))

sorted_axis(sorted, dim; rev = false) = sorted

"""
    times(x, dim=nothing)

Time coordinate of `x`: the keys of its time dimension (or of `dim`).
"""
times(x, dim = nothing) = axiskeys(x, dimnum(x, dim))

function samplingrate end

_seconds(dt::Period) = dt / Second(1)
_seconds(dt) = dt

"""
    unwrap(x)

Return the innermost object (array) of wrapped `x` with similar behavior (e.g. same size, same type, etc.)
"""
unwrap(x) = x
