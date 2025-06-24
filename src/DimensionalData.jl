function timedim(x, query=nothing)
    query = something(query, TimeDim)
    qdim = dims(x, query)
    isnothing(qdim) ? dims(x, 1) : qdim
end