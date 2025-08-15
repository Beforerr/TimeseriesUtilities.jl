module TimeseriesUtilitiesSpaceDataModelExt
import TimeseriesUtilities: times, resolution
import SpaceDataModel
using SpaceDataModel: AbstractDataVariable

times(v::AbstractDataVariable) = SpaceDataModel.times(v)
resolution(v::AbstractDataVariable) = resolution(times(v))
end
