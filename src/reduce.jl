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
tminimum(x::AbstractDimArray; query = nothing) = tminimum(times(x, query))
tmaximum(x::AbstractDimArray; query = nothing) = tmaximum(times(x, query))
targmin(x) = times(x)[argmin(x)]
targmax(x) = times(x)[argmax(x)]

"""    
    timerange(times)
    timerange(x1, xs...)

Get the time range (minimum and maximum) of time series data.

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
timerange(times) = extrema(times)

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
    # Initialize variables
    ranges = NTuple{2, eltype(times)}[]
    range_start = times[1]

    for i in 2:length(times)
        current_time = times[i]
        prev_time = times[i - 1]
        # If gap is too large, end the current range and start a new one
        if current_time - prev_time > max_dt
            push!(ranges, (range_start, prev_time))
            range_start = current_time
        end
    end
    # Add the last range
    push!(ranges, (range_start, times[end]))
    return ranges
end

"""
    find_continuous_timeranges(x, max_dt)

Find continuous time ranges for `x` (e.g. times or `DimArray`). `max_dt` is the maximum time gap between consecutive times.
"""
function find_continuous_timeranges(x, max_dt)
    isempty(x) && return []
    ts = eltype(x) <: AbstractTime ? x : times(x)
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

function time_grid(x, dt::Unitful.Quantity)
    tmin, tmax = timerange(x)
    return if dimension(dt) == Unitful.𝐓
        tmin:_2dates(dt):tmax
    elseif dimension(dt) == Unitful.𝐓^-1
        _dt = round(Nanosecond, 1 / dt)
        tmin:_dt:tmax
    else
        tmin:dt:tmax
    end
end

for (period, unit) in (
        (Dates.Week, Unitful.wk), (Dates.Day, Unitful.d), (Dates.Hour, Unitful.hr),
        (Dates.Minute, Unitful.minute), (Dates.Second, Unitful.s), (Dates.Millisecond, Unitful.ms),
        (Dates.Microsecond, Unitful.μs), (Dates.Nanosecond, Unitful.ns),
    )
    @eval _2dates(::typeof($unit)) = $period
end

_2dates(x::Unitful.Quantity) = _2dates(Unitful.unit(x))(x)
