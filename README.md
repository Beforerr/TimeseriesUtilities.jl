# TimeseriesUtilities

[![DOI](https://zenodo.org/badge/1007934281.svg)](https://doi.org/10.5281/zenodo.18634508)
[![version](https://juliahub.com/docs/General/TimeseriesUtilities/stable/version.svg)](https://juliahub.com/ui/Packages/General/TimeseriesUtilities)

Time series operations on arrays with time coordinates: selection, resampling, gap filling,
windowed statistics, derivatives, despiking, filtering and spectra.

```julia
using Pkg; Pkg.add("TimeseriesUtilities")
```

## Data model

Every operation takes a container that carries its times and returns the same kind of container:

- `Timed(data, t; dim)`: a plain array with times `t` along dimension `dim` (optional for vectors).
  It is not an `AbstractArray`, so array code outside this package throws instead of silently dropping times.
  Take it apart with `data, t = x` (or `parent(x)`, `times(x)`).
- `DimArray` (DimensionalData.jl): time is the `Ti` dimension (any `TimeDim`) or `Dim{:time}`.
- `KeyedArray` (AxisKeys.jl): time is the dimension named `:time`.

Rules shared by all functions:

- They act along the time dimension, or along the `dim` keyword: an index, a `KeyedArray` name `Symbol`, or a `DimArray` dimension.
- Times are `DateTime`/`Date` or numbers. Durations (`dt`, `every`, `window`) are `Dates.Period`s for date times, numbers for numeric times.
- Missing data is `NaN`: statistics skip it; `despike`, `tmask` and `tfill_gaps` write it.
- Time ranges `t0, t1` select the closed `[t0, t1]` or an IntervalSets interval; windows are half-open `[start, start + period)`.
- Selection and windows need sorted times (`tsort`).
- Operations that need times throw `ArgumentError` on an axis without coordinates; wrap it with `Timed`.

## Walkthrough

The data: one hour of a 3-channel signal at 1 s cadence, with two spikes and a 1-minute dropout.

```julia
using TimeseriesUtilities, Dates, Random, Statistics, LinearAlgebra

τ = 0:3599                                   # seconds
data = [sin.(2π .* τ ./ 300) cos.(2π .* τ ./ 300) τ ./ 3600] .+ 0.01 .* randn(Xoshiro(0), 3600, 3);
data[[500, 2000], 1] .= 8;                   # spikes
keep = setdiff(1:3600, 1201:1260);           # dropout 00:20:00–00:20:59
x = Timed(data[keep, :], DateTime(2024) .+ Second.(τ[keep]); dim = 1)
```

### Inspect and select

```julia
A, t = x[1:end, :]                           # indexing along time keeps the times
x[1, :]                                      # one sample: plain vector
trange(x)                                    # (2024-01-01T00:00:00, 2024-01-01T00:59:59)
cadence(t)                                   # Millisecond(1000), robust to gaps and jitter; warns about the gap
find_continuous_timeranges(t, Second(1))     # [(00:00:00, 00:19:59), (00:21:00, 00:59:59)]
tclip(x, DateTime(2024, 1, 1, 0, 10), DateTime(2024, 1, 1, 0, 20))   # copy of the samples in [t0, t1]
tselect(x, DateTime(2024, 1, 1, 0, 30))                  # nearest sample
tselect(x, DateTime(2024, 1, 1, 0, 20, 30), Second(5))   # missing: no sample within ±5 s
```

### Clean

```julia
y = despike(x, Minute(5))    # NaN where > 3 scaled MADs from the 5-minute moving median: the 2 spikes
y = tfill_gaps(y)            # insert the 60 missing times on the cadence grid, as NaN
y = tinterp_na(y)            # interpolate every NaN from its neighbors: 3600×3, no NaN
```

`outlier_mask(x, window)` returns the mask instead; `dropna(x)` drops samples containing NaN; `tmask(x, t0, t1)` blanks a known-bad interval.
Without `window`, `despike` and `outlier_mask` use the global median, which misses spikes riding on slow variations.

### Resample and align

```julia
tresample(y, Second(10))                             # linear interpolation onto a 10 s grid: 360×3
tinterp(y, DateTime(2024, 1, 1, 0, 30, 0, 500))      # one time: plain vector
v = Timed(randn(Xoshiro(1), 750, 3), DateTime(2024, 1, 1, 0, 5) .+ Second.(2:4:2998); dim = 1)  # 4 s cadence
ys, vs = tsync(y, v)          # y clipped to the common time range, v interpolated onto its times
```

Interpolating outside the data's time range throws a `DomainError` unless `extrapolation = true`.
`interp` swaps the interpolator, e.g. `tinterp(y, t_new; interp = DataInterpolations.CubicSpline)`.

### Statistics

```julia
tmean(y)                                       # per channel: plain vector
tmean(y, Minute(1))                            # 60×3, keyed by window start: [00:00, 00:01), [00:01, 00:02), …
tstat(c -> quantile(c, 0.9), y, Minute(5); period = Minute(10))  # any f, 10-minute windows every 5 minutes
tsmooth(y, Minute(1))                          # 1-minute moving mean, same size as y
tsubtract(median, y, Minute(10))               # subtract the 10-minute moving median
targmax(y)                                     # time of each channel's maximum
```

`f` receives one channel's samples in a window as a vector.
Window starts are the first time floored to `every`, then every `every` (`start_by = :datapoint` starts at the first time instead);
empty windows are kept, so `tmean` gives NaN over a dropout.

### Per-sample operations

```julia
mapsamples(norm, y)                            # |y| per time: Timed vector
mapsamples((a, b) -> cross(a, b), vs, ys)      # pairs samples by time: arguments need the same times (tsync first)
mapsamples(dot, y, Ref([0, 0, 1]))             # Ref passes a value whole to every call
tderiv(y)                                      # forward difference per second: 3599×3, keyed by the left time
```

### Frequency domain

```julia
using DSP
tfilter(y, 1 / 600, 1 / 60)                    # zero-phase band-pass between 1/600 and 1/60 Hz
pspectrum(y; nfft = 512)                       # spectrogram power: Timed 257×13×3, (frequency, time, channel)
```

### Other containers

The same calls work on `DimArray` and `KeyedArray`, and return them:

```julia
using DimensionalData
xd = DimArray(A, (Ti(t), Dim{:component}([:x, :y, :z])))
tmean(xd, Minute(30))                          # DimArray over (Ti, component)
```

Plain arrays work where no times are needed, with `dim`: statistics, `twindows`, 
`mapsamples`, `tsmooth`, `tsubtract`, `outlier_mask` and `despike`, whose `every` and `window` then count samples (windows start at the first sample)
and `targmin`/`targmax` return indices. So does a container axis without coordinates (e.g. `DimArray(A, (Ti(), Y(1:3)))`).
Other array types become containers by implementing the hooks in `src/api.jl`.

## Optional integrations

Loading a package enables:

- AxisKeys.jl, DimensionalData.jl: the containers above.
- DataInterpolations.jl, FastInterpolations.jl: `interp` constructors.
- Unitful.jl: `tresample(x, 0.1u"Hz")` or `tresample(x, 10u"s")`.
- DSP.jl: `tfilter` and `pspectrum`.

## Index

- Time coordinates: `times`, `trange`, `cadence`, `samplingrate`, `find_continuous_timeranges`, `ContinuousTimeRanges`, `tsplit`, `IntervalRange`, `tshift`
- Selection: `tclip`, `tview`, `tselect`, `tsort`, `tmask`, `tmask!`, `dropna`
- Resampling and gaps: `tinterp`, `tresample`, `tsync`, `tfill_gaps`, `tinterp_na`
- Statistics: `tstat`, `tmean`, `tmedian`, `tsum`, `tvar`, `tstd`, `tsem`, `tminimum`, `tmaximum`, `targmin`, `targmax`, `twindows`, `tsmooth`, `tsubtract`
- Cleaning: `outlier_mask`, `despike`, `despike!`
- Per sample: `mapsamples`, `tderiv`
- Frequency domain: `tfilter`, `pspectrum`

Full docstrings: `?name` in the REPL, or the [documentation](https://Beforerr.github.io/TimeseriesUtilities.jl/dev/).
