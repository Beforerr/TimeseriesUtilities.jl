# TimeseriesUtilities

[![DOI](https://zenodo.org/badge/1007934281.svg)](https://doi.org/10.5281/zenodo.18634508)
[![version](https://juliahub.com/docs/General/TimeseriesUtilities/stable/version.svg)](https://juliahub.com/ui/Packages/General/TimeseriesUtilities)

Time series operations on arrays with time coordinates: selection, resampling, gap filling,
windowed statistics, derivatives, despiking, filtering and spectra.

```julia
using Pkg; Pkg.add("TimeseriesUtilities")
```

## Data model

Every operation takes a container that carries its time coordinate and returns the same kind of container:

- `KeyedArray` (AxisKeys.jl): the time dimension is the one named `:time`.
- `DimArray` (DimensionalData.jl): the time dimension is `Ti` (any `TimeDim`) or `Dim{:time}`.
- `Timed(data, t; dim)`: wraps a plain array. Take it apart with `data, t = x` (or `parent(x)`, `times(x)`).
  It is not an `AbstractArray` intentionally, so non-package array code throws instead of silently dropping times.

Operations act along the time dimension, or along `dim` (an index, a `KeyedArray` name `Symbol`, or a `DimArray` dimension).
Times may be `DateTime`/`Date` or numbers. Periods (`dt`, `every`, `window`) are `Dates.Period`s for date times and numbers for numeric times.
Missing data is `NaN`: statistics skip NaNs, and `despike`/`tmask` write them.
Any per-sample function lifts to a container with `mapsamples`, which keeps the container and its times.

```julia
using TimeseriesUtilities, Dates, LinearAlgebra
t = DateTime(2020) .+ Second.([0, 1, 2, 5, 6])
x = Timed(rand(5, 3), t; dim = 1)                # 5 samples × 3 channels
y = tfill_gaps(x) |> tinterp_nans                 # insert samples at 3 s and 4 s, interpolate them
data, t_m = tmean(y, Second(2))                           # 2 s window means, keyed by window start
B = Timed(rand(5, 3), t; dim = 1)
v = mapsamples((e, b) -> cross(e, b) / dot(b, b), x, B)  # (x × B) / |B|² per sample: 5×3, keyed by t
```

## API

Functions are listed by call signature; keyword `dim = nothing` (time dimension) is accepted everywhere and omitted below.

### Time coordinates

- `times(x)`: time keys of `x`.
- `trange(x)`: `(tmin, tmax)`; `trange(x, ys...)`: common (intersected) range, or `nothing`.
- `cadence(t)`: robust sampling step of time vector `t` (survives gaps and jitter); `samplingrate(t)`: its inverse in Hz.
- `find_continuous_timeranges(t, max_dt)` / `ContinuousTimeRanges(t, max_dt)`: `(t0, t1)` spans split at gaps `> max_dt`.
- `tsplit(t0, t1, n_or_dt_or_PeriodType)`: split `[t0, t1)` into `n` parts, `dt`-sized parts, or calendar periods (e.g. `Month`); `IntervalRange` is its lazy form.
- `tshift(x, [t0])`: subtract `t0` (default: first time) from the times.

### Selection

- `tclip(x, t0, t1)` / `tview(x, t0, t1)`: copy / view of samples in `[t0, t1]` (times sorted; see `tsort(x; rev)`).
- `tclips(xs...)` / `tviews(xs...)`: clip / view each of `xs` to their common time range.
- `tselect(x, t, [δt])`: sample nearest to `t`; with `δt`, only within `[t-δt, t+δt]`, else `missing`.
- `tmask(x, t0, t1)`, `tmask(x, ranges)`: set samples in the range(s) to NaN (`tmask!` in place).
- `dropna(x)`: drop samples containing any NaN.
- Indexing a container along time keeps the keys: `x[2:4, :]`.

### Resampling and gaps

- `tinterp(x, t_new; interp = LinearInterpolation, extrapolation = false)`: interpolate at `t_new` (a vector, keyed by `t_new`) or at a single time (sample without the time axis).
  `interp` is any constructor `interp(u, t; kws...)` returning a callable, e.g. DataInterpolations' `CubicSpline`, or FastInterpolations' `cubic_interp` etc. (fast multi-channel path).
- `tresample(x, dt; kws...)`: `tinterp` onto a regular grid of step `dt` (also a Unitful time or frequency).
- `tsync(x, ys...)`: clip `x` to the common range and interpolate each of `ys` onto its times; returns `(x′, ys′...)`.
- `tfill_gaps(x, [dt]; fill = NaN, max_gap)`: insert missing timestamps on a `dt` (default `cadence`) grid, new samples `fill`; gaps longer than `max_gap` stay open. `tfill_gaps(t, [dt])` on a time vector returns the new times.
- `tinterp_nans(x; interp)`: interpolate NaN values from their non-NaN neighbors.

### Statistics

- `tstat(f, x)`: `f(channel_vector)` per channel over all samples; returns a scalar for a vector, else a container over the remaining axes.
- `tstat(f, x, every; period = every, start_by = :window)`: per window `[start, start + period)`, starts `every` apart from the first time floored to `every` (`:window`) or from the first time (`:datapoint`); keyed by window start. `period > every` overlaps windows. Empty windows call `f` on an empty vector.
- `tmean`, `tmedian`, `tsum`, `tvar`, `tstd`, `tsem`, `tminimum`, `tmaximum`: `tstat` with the NaN-skipping statistic, same signatures minus `f`.
- `targmin(x)`, `targmax(x)`: time of each channel's minimum/maximum, ignoring NaNs; `missing` for an all-NaN channel.
- `twindows(x, every; period, start_by)`: `(ws, starts)`, views of `x` over the same windows.
- `smooth([f], x, window)`: moving `f` (default `nanmean`) over `[c - window/2, c + window/2)` around each time (`window` may be `(before, after)`); same size as `x`.
- `tsubtract(f, x, [window])`: subtract each channel's statistic `f` (e.g. `median`, or NaNStatistics' `nanmedian` to skip NaNs), over all samples or moving over `window`.

```julia
tmean(x)                                       # one value per channel
tmaximum(x, Second(2))                         # windows [t0, t0 + 2s), starts on a regular grid
tstat(v -> count(>(0.5), v), x, Second(2); period = Second(4))  # any `f`, overlapping windows
ws, starts = twindows(x, Second(2))            # the windows themselves, as views of `x`
targmax(x)                                     # time of each channel's maximum
```

### Cleaning

- `outlier_mask(x, [window]; threshold = 3)`: `true` where a value is more than `threshold` scaled MADs from the median, globally or over a moving `window`.
- `despike(x, [window]; threshold = 3)`: set those values to NaN (`despike!` in place). Typical: `despike(x, Second(10)) |> tinterp_nans`.

### Per-sample operations

- `mapsamples(f, x, args...)`: `f(sample_of_x, samples_of_args...)` per time, stacked along time. `Ref(v)` passes `v` whole to every call.
- `proj(a, b)`, `oproj(a, b)`, `sproj(a, b)`: vector projection, rejection, scalar projection of `a` on `b` (per sample with `mapsamples`).

### Derivatives

- `tderiv(x; lazy = false)`: forward difference per second (for date times), keyed by the left time; `lazy = true` returns a lazy array.

### Frequency domain (needs `using DSP`)

- `tfilter(x, f1, [f2]; designmethod = Butterworth(2))`: zero-phase band-pass between `f1` and `f2` Hz (`f2` defaults to just below Nyquist).
- `pspectrum(x; nfft = 256, noverlap = nfft ÷ 2, window = hamming)`: spectrogram power with axes (frequency, time, remaining axes of `x`).

Full docstrings: `?name` in the REPL, or the [documentation](https://Beforerr.github.io/TimeseriesUtilities.jl/dev/).
