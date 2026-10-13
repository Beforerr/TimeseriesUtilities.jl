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

export cadence, samplingrate
export Timed
export times, targmin, targmax
export trange, find_continuous_timeranges, ContinuousTimeRanges, tsplit, IntervalRange
export tselect, tclip, tview, tmask, tmask!, tsort, tshift
export mapsamples, tsmooth, twindows
export tstat, tsum, tmean, tmedian, tstd, tsem, tvar, tmaximum, tminimum
export tderiv, tsubtract
export tinterp, tsync, tresample
export tinterp_na, tfill_gaps, dropna
export outlier_mask, despike, despike!
export tfilter, pspectrum

include("api.jl")
include("timed.jl")
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
