"""
    TimeseriesUtilities

A collection of utilities to simplify common time series analysis.
    
From data cleaning to arithmetic operations (e.g. linear algebra) to common time series operations (e.g. resampling, filtering).

## Data Cleaning

- [`find_outliers`](@ref), [`find_outliers_median`](@ref), [`find_outliers_mean`](@ref)
- [`replace_outliers`](@ref), [`replace_outliers!`](@ref)

## Query

- [`times`](@ref), [`time_grid`](@ref)
- [`timerange`](@ref), [`common_timerange`](@ref)

## (Windowed) Statistics

- Base: [`tstat`](@ref) - `tstat(f, x, [dt]; dim)`
- NaNStatistics wrappers: [`tmean`](@ref), [`tmedian`](@ref), [`tsum`](@ref), [`tvar`](@ref), [`tstd`](@ref), [`tsem`](@ref)

## Algebra

- [`tcross`](@ref), [`tdot`](@ref), [`tnorm`](@ref)
- [`tsproj`](@ref), [`tproj`](@ref), [`toproj`](@ref)
- [`tsubtract`](@ref), [`tderiv`](@ref)

## Time-Domain Operations

- [`tselect`](@ref)
- [`tclip`](@ref), [`tclips`](@ref)
- [`tview`](@ref)
- [`tmask`](@ref) and [`tmask!`](@ref)
- [`tshift`](@ref)
- [`tsplit`](@ref)
- [`tgroupby`](@ref)
- Resampling: [`tinterp`](@ref), [`tsync`](@ref)

## Time-Frequency Domain Operations

- [`tfilter`](@ref)
"""
module TimeseriesUtilities

using Base: @propagate_inbounds
using Dates
using Dates: AbstractTime
using LinearAlgebra
using StaticArrays
using NaNStatistics
using Statistics: median, median!

const SV3 = SVector{3}

export resolution, samplingrate
export times, tminimum, tmaximum, targmin, targmax
export timerange, common_timerange, time_grid, find_continuous_timeranges
export tinterp, tsync, tresample, tinterp_nans

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
export smooth, tfilter
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

"""
    tfilter(data, Wn1, Wn2=nothing; designmethod=nothing)

Bandpass filter `data` between `Wn1` and `Wn2`. The upper cutoff defaults to the Nyquist frequency.

References
- https://docs.juliadsp.org/stable/filters/
- https://www.mathworks.com/help/signal/ref/filtfilt.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html

Issues
- DSP.jl and Unitful.jl: https://github.com/JuliaDSP/DSP.jl/issues/431
"""
function tfilter end

include("compat.jl")

end
