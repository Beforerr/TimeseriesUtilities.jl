module TimeseriesUtilitiesSpaceDataModelExt
import TimeseriesUtilities: axiskeys, dimnum, tdimnum, unwrap
import SpaceDataModel as SDM
using SpaceDataModel: AbstractDataVariable

unwrap(v::AbstractDataVariable) = parent(v)
tdimnum(v::AbstractDataVariable) = SDM.tdimnum(v)
dimnum(v::AbstractDataVariable, name::Union{String, Symbol}) =
    @something findfirst(i -> SDM.name(SDM.getdim(v, i)) == name, 1:ndims(v)) throw(ArgumentError("dimension $name not found"))
axiskeys(v::AbstractDataVariable, d) = SDM.unwrap(SDM.getdim(v, d))
end
