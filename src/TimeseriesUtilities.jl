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
export trange, find_continuous_timeranges, ContinuousTimeRanges, tsplit, IntervalRange
export tselect, tclip, tclips, tview, tviews, tmask, tmask!, tsort, tshift
export mapsamples, tgroupby, smooth
export tsum, tmean, tmedian, tstd, tsem, tvar
export tderiv, tsubtract, proj, sproj, oproj
export tinterp, tsync, tresample
export tinterp_nans, tfill_gaps, dropna
export find_outliers, replace_outliers!, replace_outliers
export tfilter, pspectrum

include("api.jl")
include("utils.jl")
include("coords.jl")
include("windows.jl")
include("select.jl")
include("samples.jl")
include("stats.jl")
include("algebra.jl")
include("interp.jl")
include("gaps.jl")
include("outliers.jl")
include("spectral.jl")

end
