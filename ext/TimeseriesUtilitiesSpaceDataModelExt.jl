module TimeseriesUtilitiesSpaceDataModelExt
import TimeseriesUtilities: axiskeys, dimnum, times, unwrap
import SpaceDataModel as SDM
using SpaceDataModel: AbstractDataVariable

unwrap(v::AbstractDataVariable) = parent(v)
dimnum(v::AbstractDataVariable, dim) =
    @something dim SDM.tdimnum(v) throw(ArgumentError("$(typeof(v)) has no time dimension; pass `dim`"))
axiskeys(v::AbstractDataVariable, d) = SDM.unwrap(SDM.getdim(v, d))
times(v::AbstractDataVariable, dim = nothing) = axiskeys(v, dimnum(v, dim))
end
