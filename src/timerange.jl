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
