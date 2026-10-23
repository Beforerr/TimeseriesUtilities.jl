"""
    cadence(ts; rtol=1e-3, check=true)

Robust sampling cadence of times `ts` via IQR-trimmed modal clustering on
positive diffs, refined by the median of diffs within `rtol` of the modal center.
With `check=true`, warns on gaps or intervals off integer multiples of the cadence by more than `rtol`.
"""
cadence(ts::AbstractRange; kwargs...) = step(ts)

function cadence(ts; rtol = 1.0e-3, check = true)
    length(ts) > 1 || throw(ArgumentError("at least two time points are required"))
    T = typeof(ts[2] - ts[1])

    pos = _positive_diffs(ts)
    isempty(pos) && throw(ArgumentError("no positive time differences found"))
    sort!(pos)

    cval = _modal_cadence(pos, rtol)
    check && _check_regular(pos, cval, rtol)
    return T <: AbstractFloat ? cval : T(round(Int, cval))
end

"""
    samplingrate(ts)

Sampling rate of times `ts` in hertz: inverse [`cadence`](@ref), with numeric times taken as seconds.
"""
samplingrate(ts; kw...) = inv(_seconds(cadence(ts; kw...)))

function _positive_diffs(ts)
    @inline _float(x) = Float64(x)
    @inline _float(x::Union{Dates.Period, Dates.AbstractTime}) = Float64(Dates.value(x))

    n = length(ts) - 1
    out = Vector{Float64}(undef, n)
    any_bad = false
    @inbounds @simd for i in 1:n
        v = _float(ts[i + 1] - ts[i])
        out[i] = v
        any_bad |= (v <= 0.0)
    end
    any_bad && filter!(>(0.0), out)
    return out
end

# Pre-sorted positive diffs → fundamental cadence.
# Clustering rtol independent of precision rtol: wide enough to absorb jitter,
# narrow enough that 2Δt stays a separate cluster from Δt.
function _modal_cadence(
        pos, rtol;
        cluster_rtol = 0.05, strong_frac = 0.25, min_count = 3
    )
    n = length(pos)
    q25 = _quantile_sorted(pos, 0.25)
    q75 = _quantile_sorted(pos, 0.75)
    j = clamp(searchsortedlast(pos, q75 + 3 * (q75 - q25)), 1, n)

    p1, pj = pos[1], pos[j]
    # single-cluster check: all values within ±cluster_rtol of midpoint.
    best = if pj - p1 <= cluster_rtol * (p1 + pj)
        _median_sorted(pos, 1, j)
    else
        _smallest_strong_center(pos, j, cluster_rtol, strong_frac, min_count)
    end

    # Refine via median of sorted slice within ±rtol of best.
    i_lo = searchsortedfirst(pos, best * (1 - rtol))
    i_hi = searchsortedlast(pos, best * (1 + rtol))
    return i_lo <= i_hi ? _median_sorted(pos, i_lo, i_hi) : best
end

# Walk sorted pos[1:j], form running-mean clusters, return smallest center
# whose count clears strong threshold. Centers come out sorted, so findfirst suffices.
function _smallest_strong_center(pos, j, cluster_rtol, strong_frac, min_count)
    centers = Float64[]
    counts = Int[]
    maxk = 0
    @inbounds begin
        c = pos[1]; k = 1
        for i in 2:j
            v = pos[i]
            if v - c <= cluster_rtol * c
                k += 1
                c += (v - c) / k
            else
                push!(centers, c); push!(counts, k)
                k > maxk && (maxk = k)
                c = v; k = 1
            end
        end
        push!(centers, c); push!(counts, k)
        k > maxk && (maxk = k)
    end
    thresh = max(min_count, ceil(Int, strong_frac * maxk))
    idx = findfirst(>=(thresh), counts)
    return isnothing(idx) ? centers[1] : centers[idx]
end

function _check_regular(pos, cval, rtol)
    gap_thresh = 1.5 * cval
    bad_thresh = rtol * cval
    n_gaps = 0
    n_bad = 0
    @inbounds for v in pos
        v > gap_thresh && (n_gaps += 1)
        abs(v - round(v / cval) * cval) > bad_thresh && (n_bad += 1)
    end
    if n_gaps > 0 || n_bad > 0
        parts = String[]
        n_gaps > 0 && push!(parts, "$n_gaps gap(s)")
        n_bad > 0 && push!(parts, "$n_bad non-integer-multiple interval(s)")
        @warn "Time resolution is not approximately constant ($(join(parts, ", ")))"
    end
    return nothing
end

"""
    targmin(x; dim=nothing)
    targmax(x; dim=nothing)

Time of the minimum (maximum) of each channel of `x`, ignoring NaNs; `missing` for an all-NaN channel.
"""
targmin(x; dim = nothing) = _targ(nanargmin, x, dim)
targmax(x; dim = nothing) = _targ(nanargmax, x, dim)

function _targ(argf::F, x, dim) where {F}
    t = _keys_or_indices(x, resolve_dim(x, dim))
    # `nanargmax` returns 1 for an empty or all-NaN slice
    return tstat(x; dim) do v
        i = isempty(v) ? 0 : argf(v)
        i == 0 || isnan(v[i]) ? missing : t[i]
    end
end

"""
    trange(x; dim=nothing)
    trange(x1, xs...; dim=nothing)

Time range `(tmin, tmax)` of `x`. With several arguments, their common time range
(intersection), or `nothing` if they do not overlap.
"""
trange(x; dim = nothing) = extrema(times(x, dim))

function trange(x1, xs...; dim = nothing)
    t0, t1 = trange(x1; dim)
    for x in xs
        _t0, _t1 = trange(x; dim)
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

"""
    ContinuousTimeRanges(times, max_dt)

Iterator of `(t_start, t_stop)` spans over `times`, splitting wherever consecutive values differ by more than `max_dt`.
"""
struct ContinuousTimeRanges{T, P}
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
struct IntervalRange{T, D}
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

# Handle numeric offsets for datetime-like types (default to Unix epoch)
struct TimeOffsets{T, A <: AbstractArray} <: AbstractVector{T}
    offsets::A
    t0::T
end

TimeOffsets(offsets) = TimeOffsets(offsets, Dates.unix2datetime(0))

Base.size(to::TimeOffsets) = size(to.offsets)
function Base.getindex(to::TimeOffsets, i::Int)
    _add(t0, dt) = t0 + dt
    _add(t0::Dates.AbstractTime, dt::Number) = t0 + round(Nanosecond(round(Int, 1.0e9 * dt)), typeof(t0 - t0))
    return _add(to.t0, to.offsets[i])
end

"""
    tsplit(t0, t1, n::Int)
    tsplit(t0, t1, dt)
    tsplit(t0, t1, dtType::Type{<:Period})

Split the range from `t0` to `t1` into `n` parts, `dt`-sized parts, or by period type (e.g., Month, Day).
"""
tsplit(t0, t1, dt) = collect(IntervalRange(t0, t1, dt))
tsplit((t0, t1), arg) = tsplit(t0, t1, arg)
