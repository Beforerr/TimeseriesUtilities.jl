# TimeseriesUtilities

[![DOI](https://zenodo.org/badge/1007934281.svg)](https://doi.org/10.5281/zenodo.18634508)
[![version](https://juliahub.com/docs/General/TimeseriesUtilities/stable/version.svg)](https://juliahub.com/ui/Packages/General/TimeseriesUtilities)

[![Build Status](https://github.com/Beforerr/TimeseriesUtilities.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Beforerr/TimeseriesUtilities.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Beforerr/TimeseriesUtilities.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Beforerr/TimeseriesUtilities.jl)
[![](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


A collection of utilities to simplify common time series analysis.

**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("TimeseriesUtilities")`

## Quickstart

Use `tinterp` directly with arrays when time coordinates are passed separately:

```julia
using TimeseriesUtilities

t = [0.0, 1.0, 2.0, 3.0]
x = [0.0, 1.0, 4.0, 9.0]

tinterp(x, t, [0.5, 1.5, 2.5])
# 3-element Vector{Float64}: [0.5, 2.5, 6.5]
```

For time-aware arrays, use a container package such as DimensionalData.jl. The
time axis is inferred, and operations rebuild the same kind of container:

```julia
using Dates
using DimensionalData
using TimeseriesUtilities

t = DateTime(2020, 1, 1):Hour(1):DateTime(2020, 1, 1, 5)
x = DimArray([0.0, 1.0, 4.0, 9.0, 16.0, 25.0], (Ti(t),))

tclip(x, DateTime(2020, 1, 1, 1), DateTime(2020, 1, 1, 3))
# 3-element DimArray with Ti axis at 01:00, 02:00, 03:00

tinterp(x, DateTime(2020, 1, 1, 2, 30))
# 6.5

tresample(x, Minute(30))
# DimArray resampled onto a regular 30-minute time grid
```

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg?logo=julia)](https://Beforerr.github.io/TimeseriesUtilities.jl/dev/)
