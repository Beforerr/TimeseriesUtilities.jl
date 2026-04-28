# A function-like wrapper for time series interpolation that handles time type conversions.
# a workaround for `Time` type: https://github.com/SciML/DataInterpolations.jl/issues/436
struct Tinterp{T, F}
    interp::F
    tType::Type{T}
end

Tinterp(u, t, interp; kws...) = Tinterp(interp(u, rawview(t); kws...), eltype(t))
Tinterp(interp) = (u, t; kws...) -> Tinterp(u, t, interp; kws...)
(ti::Tinterp{T})(t) where {T} = ti.interp(rawview(T(t)))
