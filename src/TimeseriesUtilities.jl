module TimeseriesUtilities

@doc let path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end TimeseriesUtilities

using Base: @propagate_inbounds
using Dates
using Dates: AbstractTime
using LinearAlgebra
using NaNStatistics
using Statistics: mean, median, median!, quantile, middle

export cadence, resolution, samplingrate
export times, tminimum, tmaximum, targmin, targmax
export timerange, common_timerange, time_grid, find_continuous_timeranges
export tinterp, tsync, tresample, tinterp_nans, tfill_gaps

# Time operations
export tselect, tclip, tclips, tview, tviews, tmask, tmask!, tsort, tshift
# Linear Algebra
export proj, sproj, oproj
export tdot, tcross, tnorm, tproj, tsproj, toproj, tnorm_combine
export tgroupby
# Statistics
export tsum, tmean, tmedian, tstd, tsem, tvar
# Derivatives
export tderiv, tsubtract

# Data cleaning
export smooth
export dropna
export find_outliers, replace_outliers!, replace_outliers

export tsplit, IntervalRange

include("api.jl")
include("sliding.jl")
include("timerange.jl"); export ContinuousTimeRanges
include("timeseries.jl")
include("operations.jl")
include("groupby.jl")
include("reduce.jl")
include("stats.jl")
include("algebra.jl")
include("lazyoperations.jl")
include("interp.jl")
include("outliers.jl")
include("utils.jl")
include("timefrequency.jl"); export tfilter, pspectrum

include("compat.jl")

end
