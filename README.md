# TimeseriesUtilities

[![DOI](https://zenodo.org/badge/1007934281.svg)](https://doi.org/10.5281/zenodo.18634508)
[![version](https://juliahub.com/docs/General/TimeseriesUtilities/stable/version.svg)](https://juliahub.com/ui/Packages/General/TimeseriesUtilities)

A collection of utilities to simplify common time series analysis:
from data cleaning to arithmetic operations (e.g. linear algebra) to common time series operations (e.g. resampling, filtering).

**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("TimeseriesUtilities")`

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg?logo=julia)](https://Beforerr.github.io/TimeseriesUtilities.jl/dev/)

Most of the utilities operate on the time dimension by default, but you can specify other dimensions using the `dim` or `query` parameter.

## Data Cleaning

- `find_outliers`, `find_outliers_median`, `find_outliers_mean`
- `replace_outliers`, `replace_outliers!`

## Query

- `times`
- `trange`
- `cadence` — robust modal-cluster Δt, survives gaps, jitter, and dropout

## (Windowed) Statistics

- Base: `tstat` - `tstat(f, x, [dt]; dim)`
- NaNStatistics wrappers: `tmean`, `tmedian`, `tsum`, `tvar`, `tstd`, `tsem`

## Algebra

- `tcross`, `tdot`, `tnorm`
- `tsproj`, `tproj`, `toproj`
- `tsubtract`, `tderiv`

## Per-Sample Operations

- `eachsample(A, [ts]; dim)` — iterate `(x, t)` over slices along time. `ts` is either one time per sample or a single value shared by all samples. It defaults to the time keys of `A`; plain arrays have none and must pass it.
- `mapsamples(f, A, [ts]; dim)` — apply `f(x, t)` per sample and stack results along time (at position `dim`, or last when the result has fewer axes)

## Time-Domain Operations

- `tselect`
- `tclip`, `tclips`
- `tview`
- `tmask` and `tmask!`
- `tshift`
- `tsplit`
- `tgroupby`
- Resampling: `tinterp`, `tsync`
- `tfill_gaps` — insert missing timestamps as `fill` (NaN), originals untouched

## Time-Frequency Domain Operations

- `tfilter`, `pspectrum`
