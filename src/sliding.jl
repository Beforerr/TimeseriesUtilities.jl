struct SlidingWindow{T, D, C, B, A} <: AbstractVector{T}
    data::D
    coords::C
    before::B
    after::A
end

function SlidingWindow(data::D, coords::C, before::B, after::A) where {D, C, B, A}
    T = Core.Compiler.return_type(view, Tuple{D, UnitRange{Int}})
    return SlidingWindow{T, D, C, B, A}(data, coords, before, after)
end

Base.size(sw::SlidingWindow) = (length(sw.data),)

Base.@propagate_inbounds function Base.getindex(sw::SlidingWindow, i::Int)
    coord = sw.coords[i]
    lo = searchsortedfirst(sw.coords, coord - sw.before)
    hi = searchsortedfirst(sw.coords, coord + sw.after) - 1
    return view(sw.data, lo:hi)
end