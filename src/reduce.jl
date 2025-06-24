"""
    tminimum(x)

Get the minimum timestamp of `x`.
"""
function tminimum end

"""
    tmaximum(x)

Get the maximum timestamp of `x`.
"""
function tmaximum end

tminimum(x) = minimum(x)
tmaximum(x) = maximum(x)
tminimum(x::AbstractDimArray; query=nothing) = tminimum(times(x, query))
tmaximum(x::AbstractDimArray; query=nothing) = tmaximum(times(x, query))

timerange(times) = _extrema(times)

_extrema(x) = extrema(x)
function _extrema(x::Array{T}) where {T <: Union{Date, DateTime, Int}}
    return reinterpret.(T, vextrema(reinterpret(Int, x)))
end

timerange(x1, xs...) = common_timerange(x1, xs...)

"""
    common_timerange(arrays)

Get the common time range (intersection) across multiple arrays.
If there is no overlap, returns nothing.
"""
function common_timerange(x1, xs...)
    t0, t1= timerange(x1)
    for x in xs
        _t0, _t1 = timerange(x)
        t0 = max(t0, _t0)
        t1 = min(t1, _t1)
        t0 > t1 && return nothing
    end
    return t0, t1
end