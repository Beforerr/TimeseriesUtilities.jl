abstract type AbstractTimeRanges end

"""
    ContinuousTimeRanges(times, max_dt)

Iterator of `(t_start, t_stop)` spans over `times`, splitting wherever consecutive values differ by more than `max_dt`.
"""
struct ContinuousTimeRanges{T, P} <: AbstractTimeRanges
    times::T
    max_dt::P
end

Base.IteratorSize(::Type{<:ContinuousTimeRanges}) = Base.SizeUnknown()
Base.eltype(::Type{<:ContinuousTimeRanges{T}}) where {T} = NTuple{2, eltype(T)}

function Base.iterate(iter::ContinuousTimeRanges, start_idx = 1)
    times = iter.times
    start_idx > length(times) && return nothing

    range_start = times[start_idx]
    prev_time = range_start

    for i in (start_idx + 1):length(times)
        current_time = times[i]
        if current_time - prev_time > iter.max_dt
            return (range_start, prev_time), i
        end
        prev_time = current_time
    end
    return (range_start, last(times)), length(times) + 1
end


# References:
# - pandas.interval_range : https://pandas.pydata.org/docs/reference/api/pandas.interval_range.html
# - Arrow span_range : https://arrow.readthedocs.io/en/latest/api-guide.html
"""
    IntervalRange{T, D}

Lazy iterator of `(t_start, t_end)` pairs splitting `[t0, t1)` into `dt`-sized windows.
"""
struct IntervalRange{T, D} <: AbstractTimeRanges
    t0::T
    t1::T
    dt::D
end

IntervalRange(t0::T1, t1::T2, dtType::Type{<:Period}) where {T1, T2} =
    IntervalRange(t0, t1, dtType(1))

function IntervalRange(t0::T1, t1::T2, dt::D) where {T1, T2, D}
    T = promote_type(T1, T2)
    return IntervalRange(T(t0), T(t1), dt)
end
function IntervalRange(t0::T1, t1::T2, n::Int) where {T1, T2}
    @assert n > 0
    return IntervalRange(t0, t1, (t1 - t0) /ₜ n)
end


Base.eltype(::Type{<:IntervalRange{T}}) where {T} = NTuple{2, T}

const CalendarPeriod = Union{Dates.Month, Dates.Quarter, Dates.Year}

_months(dt::Dates.Month) = Dates.value(dt)
_months(dt::Dates.Quarter) = 3Dates.value(dt)
_months(dt::Dates.Year) = 12Dates.value(dt)

Base.length(s::IntervalRange) = s.t0 >= s.t1 ? 0 : ceil(Int, (s.t1 - s.t0) / s.dt)

# Month/Quarter/Year: estimate from calendar month ordinals, then correct for day/time.
function Base.length(s::IntervalRange{T, D}) where {T, D <: CalendarPeriod}
    s.t0 >= s.t1 && return 0
    step = _months(s.dt)
    step <= 0 && throw(ArgumentError("dt must be positive"))
    n = cld(12year(s.t1) + month(s.t1) - (12year(s.t0) + month(s.t0)), step)
    return s.t0 + n * s.dt < s.t1 ? n + 1 : n
end

function Base.iterate(s::IntervalRange, state = (1, length(s)))
    i, n = state
    i > n && return nothing
    return s[i], (i + 1, n)
end

function Base.getindex(s::IntervalRange{T}, i::Int) where {T}
    @boundscheck (1 <= i <= length(s)) || throw(BoundsError(s, i))
    t_start = s.t0 + (i - 1) * s.dt
    return (t_start, min(s.t0 + i * s.dt, s.t1))
end
