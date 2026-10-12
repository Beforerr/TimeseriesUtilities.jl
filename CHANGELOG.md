# Changelog

## Unreleased

### Added

- `ws, starts = twindows(x, every; period, start_by)`: views of `x` over the windows of windowed `tstat`.
- `Timed(data, t; dim)` attaches times to a plain array; every operation accepts and returns it. Take results apart with `data, t = x`.

### Changed

- **Breaking**: Explicit-times forms are removed; wrap plain arrays with `Timed` instead: `tinterp(A, t, t_new)`, `tresample(A, t, dt)`, `tfill_gaps(A, t, [dt])`, `smooth(A, t, window)`, `tderiv(A, t)`, `tfilter(A, t, Wn1, Wn2)`, `pspectrum(A, t)`.
- **Breaking**: Operations on plain arrays that need time keys throw `ArgumentError` instead of `MethodError`.
- **Breaking**: Containers that unwrap to other data but lack a `rebuild` method throw instead of returning plain data.
- **Breaking**: `smooth(x, window)` windows over the axis keys of every container (it used sample indices for `KeyedArray`); plain arrays still use sample indices.
- **Breaking**: `tinterp(A, t)` at a single time on a multidimensional `A` returns the sample without the time axis.
- **Breaking**: Windowed statistics take `tstat(f, x, every; period = every, start_by = :window)` (same for `tmean`, `tsum`, ...); `period > every` overlaps windows. Empty windows are kept: starts form a regular grid and `f` is called on an empty vector.
- **Breaking**: `tmaximum(x, [every])` and `tminimum` are statistics of the data like `tmean`.
- **Breaking**: `targmin(x)` and `targmax(x)` return one time per channel.
- **Breaking**: `tstat(f, x)` calls `f` on each channel's samples as a vector instead of `f(data; dim)`. `tstat` is exported.
- **Breaking**: Outlier functions are replaced by `despike(A, [window]; threshold = 3)` (and `despike!`), which sets outliers to NaN.
- Full reductions (`tmean(x)` etc.) keep the remaining axes of `KeyedArray` as they do for `DimArray`.
- **Breaking**: `cadence`, `resolution`, `samplingrate`, and `find_continuous_timeranges` take time coordinates only now; pass `times(x)` for a container.
- **Breaking**: Rename `timerange(x)` to `trange(x)`; it accepts containers only. `common_timerange(xs...)` becomes `trange(xs...)`.
- **Breaking**: `times(x)` on a plain array errors instead of returning `x`; so do `targmin` and `targmax`.
- **Breaking**: Plain arrays with more than one dimension have no default time dimension: calls without `dim` throw `ArgumentError`. Plain vectors are unaffected.
- **Breaking**: `mapsamples(f, A, args...)` calls `f` with the matching sample of each argument, each sliced along its own time dimension. Time is no longer passed implicitly: `mapsamples(f, A)` becomes `mapsamples(f, A, times(A))`. Keyed arguments must share the time keys of `A`.
- **Breaking**: Remove the `trange` keyword of `tclips` and `tviews`; use `tclip.(xs, Ref(tr))`.
- **Breaking**: `samplingrate(ts)` is defined without Unitful and returns hertz as a plain number.
- `times(x, dim)` is generic over `axiskeys` and `tdimnum`; container backends no longer define `times`.
- **Breaking**: `smooth(f, x, window)` takes the moving statistic first, called on each channel's window as a vector; the `op` keyword is removed. `smooth(x, window)` keeps the `nanmean` default.
- **Breaking**: `tsubtract(f, x, [window])` requires the statistic `f`, called on each channel's samples as a vector instead of `f(data; dims)`; with `window` it subtracts the moving statistic of `smooth(f, x, window)`.
- **Breaking**: Rename `tinterp_nans` to `tinterp_na`, matching `dropna`.

### Fixed

- `tderiv` along a dimension other than the first no longer divides by time differences laid along the first.

### Removed

- **Breaking**: The SpaceDataModel extension; SpaceDataModel is to implement the container interface for its own variables.
- **Breaking**: `tgroupby`, `groupby_dynamic` and `GroupByDynamic`; use `tstat(f, x, every; period)` or `twindows`.
- **Breaking**: `resolution`; use `cadence`.
- **Breaking**: `tview(f::Function, t0, t1)`.
- **Breaking**: `Tinterp`; external interpolators passed as `interp` receive raw time values (`Dates.value`) directly.
- **Breaking**: `time_grid`; use `tresample(A, dt)` or `range(trange(x)...; step = dt)`.
- **Breaking**: `tnorm`, `tdot`, `tcross`, `tsproj`, `tproj`, `toproj`, `tnorm_combine`; use `mapsamples(norm, B)`, `mapsamples(cross, E, B)`, `mapsamples(oproj, A, B)`, ...
- **Breaking**: `tfill_gaps(t, [dt])` on time vectors; use `tfill_gaps(x, [dt])` on a container.

## v0.2.0

### Added

- Move keyed time operations to the generic axis API, enabling `AxisKeys.KeyedArray` support for operations such as `tsort`, `tclip`, `tview`, `tselect`, `tmask`, `tshift`, and grouped statistics.

### Changed

- **Breaking**: Move DimensionalData and Unitful integrations from hard dependencies to package extensions.
- **Breaking**: `tderiv` now returns plain floating-point derivatives per second for floating-point data with `Date`/`DateTime` coordinates, so Unitful is no longer required.
- **Breaking**: Remove container-target interpolation such as `tinterp(A, B::AbstractDimArray)`; pass target coordinates explicitly with `tinterp(A, times(B))`.
