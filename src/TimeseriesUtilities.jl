module TimeseriesUtilities

using Dates
using DimensionalData
using DimensionalData.Dimensions
using DimensionalData.Lookups
using VectorizedStatistics, NaNStatistics
using Unitful

export timerange, common_timerange

# Time operations
export tselect, tclip, tclips, tview, tviews, tmask, tmask!, tsort, tshift
# Linear Algebra
export proj, sproj, oproj
export tdot, tcross, tnorm, tproj, tsproj, toproj
export tgroupby
# Statistics
export tsum, tmean, tmedian, tstd, tsem, tvar
# Derivatives
export tderiv, tsubtract

export find_outliers, replace_outliers!, replace_outliers

include("timeseries.jl")
include("operations.jl")
include("groupby.jl")
include("reduce.jl")
include("stats.jl")
include("methods.jl")
include("lazyoperations.jl")
include("outliers.jl")
include("utils.jl")
include("DimensionalData.jl")

end
