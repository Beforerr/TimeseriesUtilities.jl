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
    timerange(x)
    timerange(x1, xs...)

Get the time range of time series data `x`.

For a single argument, returns a tuple `(tmin, tmax)` containing the minimum and maximum times.
For multiple arguments, returns the common time range (intersection) across all arrays - equivalent to `common_timerange(x1, xs...)`.

# Examples
```julia
# Single time series
times = [1, 2, 3, 4, 5]
timerange(times)  # (1, 5)

# Multiple time series - find common range
x1_times = [1, 2, 3, 4]
x2_times = [2, 3, 4, 5]
timerange(x1_times, x2_times)  # (2, 4)
```

See also: [`common_timerange`](@ref), [`tminimum`](@ref), [`tmaximum`](@ref)
"""
timerange(x) = extrema(times(x))

_median(x) = median(x)

for f in (:median, :median!)
    _f = Symbol(:_, f)
    @eval $_f(x::AbstractArray{T}) where {T <: Union{Date, DateTime, Period, Int}} = T(round(Int, $f(reinterpret(Int, x))))
end

# function _median(x::AbstractArray{T}) where {T<:Union{Date,DateTime,Period}}
#     return T(round(Int, median(reinterpret(Int, x))))
# end

timerange(x1, xs...) = common_timerange(x1, xs...)

"""
    common_timerange(x1, xs...)

Get the common time range (intersection) across multiple arrays.
If there is no overlap, returns nothing.
"""
function common_timerange(x1, xs...)
    t0, t1 = timerange(x1)
    for x in xs
        _t0, _t1 = timerange(x)
        t0 = max(t0, _t0)
        t1 = min(t1, _t1)
        t0 > t1 && return nothing
    end
    return t0, t1
end

function _find_continuous_timeranges(times, max_dt)
    T = eltype(times)
    ranges = NTuple{2, T}[]
    isempty(times) && return ranges

    prev_time = range_start = first(times)

    for current_time in @view times[2:end]
        if current_time - prev_time > max_dt
            push!(ranges, (range_start, prev_time))
            range_start = current_time
        end
        prev_time = current_time
    end
    push!(ranges, (range_start, last(times)))
    return ranges
end

"""
    find_continuous_timeranges(x, max_dt)

Find continuous time ranges for `x`, where `max_dt` is the maximum time gap between consecutive times.
"""
function find_continuous_timeranges(x, max_dt)
    ts = times(x)
    return issorted(ts) ?
           _find_continuous_timeranges(ts, max_dt) :
           _find_continuous_timeranges(sort(ts), max_dt)
end


"""
    time_grid(x, dt)

Create a time grid from the minimum to maximum time in `x` with the step size `dt`.

# Examples
```julia
# Create hourly time grid
time_grid(x, Hour(1))
time_grid(x, 1u"hr")

# Create 1-s intervals
time_grid(x, Second(1))
time_grid(x, 1u"second")
time_grid(x, 1u"Hz")
```
"""
function time_grid(x, dt)
    tmin, tmax = timerange(x)
    return tmin:dt:tmax
end
