"""
    TimeseriesUtilities

A collection of utilities to simplify common time series analysis.
    
From data cleaning to arithmetic operations (e.g. linear algebra) to common time series operations (e.g. resampling, filtering).

## Data Cleaning

- [`find_outliers`](@ref), [`find_outliers_median`](@ref), [`find_outliers_mean`](@ref)
- [`replace_outliers`](@ref), [`replace_outliers!`](@ref)

## (Windowed) Statistics

- [`tstat`](@ref)
- [`tmean`](@ref)
- [`tmedian`](@ref)
- [`tsum`](@ref)
- [`tvar`](@ref)
- [`tstd`](@ref)
- [`tsem`](@ref)

## Algebra

- [`tcross`](@ref)
- [`tdot`](@ref)
- [`tnorm`](@ref)
- [`tsproj`](@ref), [`tproj`](@ref), [`toproj`](@ref)
- [`tsubtract`](@ref)
- [`tderiv`](@ref)

## Time-Domain Operations

- [`tselect`](@ref)
- [`tclip`](@ref), [`tclips`](@ref)
- [`tview`](@ref)
- [`tmask`](@ref) and [`tmask!`](@ref)
- [`tshift`](@ref)
- [`tsplit`](@ref)
- [`tgroupby`](@ref)

## Time-Frequency Domain Operations

- [`tfilter`](@ref)
"""
module TimeseriesUtilities

using Dates
using Dates: AbstractTime
using DimensionalData
using DimensionalData.Dimensions
using DimensionalData.Lookups
using LinearAlgebra
using LinearAlgebra: norm2
using StaticArrays
using VectorizedStatistics, NaNStatistics
using Unitful
using DataInterpolations
using DSP
using RollingWindowArrays

const SV3 = SVector{3}

export resolution, samplingrate
export tminimum, tmaximum, timerange, common_timerange

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
export dropna, rectify
export find_outliers, replace_outliers!, replace_outliers

export tsplit

include("timeseries.jl")
include("operations.jl")
include("groupby.jl")
include("reduce.jl")
include("stats.jl")
include("algebra.jl")
include("lazyoperations.jl")
include("outliers.jl")
include("utils.jl")
include("DimensionalData.jl")

end
