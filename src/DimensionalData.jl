dimnum(x, query) = DimensionalData.dimnum(x, something(query, TimeDim))
times(x::AbstractDimArray, args...) = parent(lookup(timedim(x, args...)))

function timedim(x, query=nothing)
    query = something(query, TimeDim)
    qdim = dims(x, query)
    isnothing(qdim) ? dims(x, 1) : qdim
end

timerange(times::DimensionalData.Sampled) = timerange(parent(times))
timerange(times::Dimension) = timerange(parent(times))
timerange(x::AbstractDimArray) = timerange(times(x))

function groupby_dynamic(x::Dimension, args...; kwargs...)
    return groupby_dynamic(parent(lookup(x)), args...; kwargs...)
end