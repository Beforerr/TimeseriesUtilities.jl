module TimeseriesUtilitiesSpaceDataModelExt
import TimeseriesUtilities: times
import SpaceDataModel
using SpaceDataModel: AbstractDataVariable

times(v::AbstractDataVariable) = SpaceDataModel.times(v)
end
