module TimeseriesUtilitiesStaticArraysExt

using TimeseriesUtilities
using StaticArrays
import TimeseriesUtilities: LazySlices

# Specialize LazySlices for numeric arrays to use SArray slices.
# SArray arithmetic returns SArray, satisfying the round-trip requirement of splines.
function TimeseriesUtilities.LazySlices(A::AbstractArray{Tv, N}, d::Int) where {Tv <: Number, N}
    other = ntuple(i -> size(A, i < d ? i : i + 1), N - 1)
    S = SArray{Tuple{other...}, Tv, N - 1, prod(other)}
    return TimeseriesUtilities.LazySlices{S, typeof(A), d}(A)
end

@inline function Base.getindex(s::LazySlices{S, A, d}, k::Int) where {S <: SArray, A, d}
    slice_sz = size(S)
    return S(
        ntuple(Val(length(S))) do p
            ci = CartesianIndices(slice_sz)[p]
            full_idx = ntuple(Val(ndims(A))) do j
                j < d ? ci[j] : j == d ? k : ci[j - 1]
            end
            @inbounds s.data[full_idx...]
        end
    )
end

end
