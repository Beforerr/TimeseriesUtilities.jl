function timedim(x, query=nothing)
    query = something(query, TimeDim)
    qdim = dims(x, query)
    isnothing(qdim) ? dims(x, 1) : qdim
end

timerange(times::DimensionalData.Sampled) = timerange(parent(times))
timerange(times::Dimension) = timerange(parent(times))
timerange(x::AbstractDimArray) = timerange(times(x))