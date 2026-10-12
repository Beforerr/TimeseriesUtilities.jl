module TimeseriesUtilitiesDataInterpolationsExt

using DataInterpolations: AbstractInterpolation, ExtrapolationType
import TimeseriesUtilities as TU

TU._interp_kws(::Type{<:AbstractInterpolation}; extrapolation = nothing, kws...) =
    isnothing(extrapolation) ? kws : (; extrapolation = _extrap(extrapolation), kws...)

_extrap(e::Bool) = e ? ExtrapolationType.Extension : ExtrapolationType.None
_extrap(e) = e

end
