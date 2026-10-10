"""
    tminimum(x)

Get the minimum timestamp of `x`.
"""
tminimum(x) = minimum(times(x))

"""
    tmaximum(x)

Get the maximum timestamp of `x`.
"""
tmaximum(x) = maximum(times(x))

targmin(x) = times(x)[argmin(x)]
targmax(x) = times(x)[argmax(x)]

"""
    trange(x)
    trange(x1, xs...)

Time range `(tmin, tmax)` of `x`. With several arguments, their common time range
(intersection), or `nothing` if they do not overlap.
"""
trange(x) = extrema(times(x))

function trange(x1, xs...)
    t0, t1 = trange(x1)
    for x in xs
        _t0, _t1 = trange(x)
        t0 = max(t0, _t0)
        t1 = min(t1, _t1)
        t0 > t1 && return nothing
    end
    return t0, t1
end

"""
    find_continuous_timeranges(ts, max_dt)

Find continuous time ranges in times `ts`, where `max_dt` is the maximum time gap between consecutive times.
"""
function find_continuous_timeranges(ts, max_dt)
    return collect(
        issorted(ts) ?
            ContinuousTimeRanges(ts, max_dt) :
            ContinuousTimeRanges(sort(ts), max_dt)
    )
end

function _time_grid(ts, dt)
    tmin, tmax = extrema(ts)
    return tmin:dt:tmax
end
