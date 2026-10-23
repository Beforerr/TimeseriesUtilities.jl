# Changelog

## v0.4.0

### Changed

- **Breaking**: a dimension without coordinates has `axiskeys` `nothing` (plain arrays included); `times` throws for one.
- **Breaking**: backends' `rebuild` receives `nothing` keys for a dimension without coordinates, e.g. those `mapsamples` adds for each sample's result; never for a dimension of `x` that has them.
- **Breaking**: `mapsamples` throws `DimensionMismatch` when an argument whose samples have the shape of those of `A` has different coordinates on a sample axis.

## v0.3.0

Every operation now takes a container that carries its times (`KeyedArray`, `DimArray`, or the new `Timed`) and returns the same kind of container.
Explicit-times forms, per-sample vector algebra and one-off helpers are removed in favor of a smaller set of composable operations.

### Added

- `Timed(data, t; dim)` attaches times to a plain array; every operation accepts and returns it. Take results apart with `data, t = x`.
- `tstat(f, x, [every])` is exported: any per-channel statistic, over all samples or per window.
- `ws, starts = twindows(x, every; period, start_by)`: views of `x` over the windows of windowed `tstat`.
- `outlier_mask(A, [window]; threshold)`, `despike(A, [window]; threshold)` and `despike!`: MAD outliers, and setting them to NaN.

### Changed

#### Data model

- **Breaking**: Explicit-times forms are removed; wrap plain arrays with `Timed` instead: `tinterp(A, t, t_new)`, `tresample(A, t, dt)`, `tfill_gaps(A, t, [dt])`, `smooth(A, t, window)`, `tderiv(A, t)`, `tfilter(A, t, Wn1, Wn2)`, `pspectrum(A, t)`.
- **Breaking**: Operations on plain arrays that need time keys throw `ArgumentError` instead of `MethodError`.
- **Breaking**: Containers that unwrap to other data but lack a `rebuild` method throw instead of returning plain data.
- **Breaking**: Plain arrays with more than one dimension have no default time dimension: calls without `dim` throw `ArgumentError`. Plain vectors are unaffected.
- **Breaking**: `times(x)` on a plain array errors instead of returning `x`.
- **Breaking**: `dimnum(x, dim)` errors on `nothing`; operations resolve their `dim` keyword with `resolve_dim(x, dim)`, where `nothing` selects the time dimension.
- **Breaking**: `tdimnum` of a plain vector is `nothing`, so it means a declared time dimension; operations on plain vectors still default to their only dimension. Containers that are vectors need a `tdimnum` method to be used without `dim`.
- `times(x, dim)` is generic over `axiskeys` and `tdimnum`; container backends no longer define `times`.

#### Time coordinates

- **Breaking**: `cadence`, `samplingrate`, and `find_continuous_timeranges` take time coordinates only; pass `times(x)` for a container.
- **Breaking**: `samplingrate(ts)` is defined without Unitful and returns hertz as a plain number.
- **Breaking**: Rename `timerange(x)` to `trange(x)`; it accepts containers only. `common_timerange(xs...)` becomes `trange(xs...)`.


#### Statistics and windows

- **Breaking**: `tstat(f, x)` calls `f` on each channel's samples as a vector instead of `f(data; dim)`.
- **Breaking**: Windowed statistics take `tstat(f, x, every; period = every, start_by = :window)` (same for `tmean`, `tsum`, ...); `period > every` overlaps windows. Empty windows are kept: starts form a regular grid and `f` is called on an empty vector.
- **Breaking**: `tmaximum(x, [every])` and `tminimum` are statistics of the data like `tmean`.
- **Breaking**: `targmin(x)` and `targmax(x)` return one time per channel, ignoring NaNs.
- **Breaking**: Rename `smooth` to `tsmooth`, matching the other time operations. `tsmooth(f, x, window)` takes the moving statistic first, called on each channel's window as a vector; the `op` keyword is removed. `tsmooth(x, window)` keeps the `nanmean` default.
- **Breaking**: `tsubtract(f, x, [window])` requires the statistic `f`, called on each channel's samples as a vector instead of `f(data; dims)`; with `window` it subtracts the moving statistic of `tsmooth(f, x, window)`.

#### Per-sample operations

- **Breaking**: `mapsamples(f, A, args...)` calls `f` with the matching sample of each argument, each sliced along its own time dimension. Time is no longer passed implicitly: `mapsamples(f, A)` becomes `mapsamples(f, A, times(A))`. Keyed arguments must share the time keys of `A`.
- **Breaking**: `proj`, `sproj` and `oproj` are no longer exported.

### Fixed

- `tderiv` along a dimension other than the first no longer divides by time differences laid along the first.

### Removed / Renamed

- **Breaking**: Rename `tinterp_nans` to `tinterp_na`, matching `dropna`.
- **Breaking**: `find_outliers`, `replace_outliers` and `replace_outliers!`; use `outlier_mask` and `despike`, then `tinterp_na` to fill. The `:mean` detection method is dropped, and `outlier_mask` is global unless given a `window`.
- **Breaking**: `tgroupby`, `groupby_dynamic` and `GroupByDynamic`; use `tstat(f, x, every; period)` or `twindows`.
- **Breaking**: `eachsample`; use `mapsamples`, or `zip(eachslice(parent(x); dims), times(x))`.
- **Breaking**: `tnorm`, `tdot`, `tcross`, `tsproj`, `tproj`, `toproj`, `tnorm_combine`; use `mapsamples(norm, B)`, `mapsamples(cross, E, B)`, `mapsamples(TimeseriesUtilities.oproj, A, B)`, ...
- **Breaking**: `resolution`; use `cadence`.
- **Breaking**: `time_grid`; use `tresample(A, dt)` or `range(trange(x)...; step = dt)`.
- **Breaking**: `tfill_gaps(t, [dt])` on time vectors; use `tfill_gaps(x, [dt])` on a container.
- **Breaking**: `tclips` and `tviews`; use `tclip.(xs, Ref(trange(xs...)))` or `tview.(...)`.
- **Breaking**: `tview(f::Function, t0, t1)`.
- **Breaking**: `Tinterp`; external interpolators passed as `interp` receive raw time values (`Dates.value`) directly.

## v0.2.0

### Added

- Move keyed time operations to the generic axis API, enabling `AxisKeys.KeyedArray` support for operations such as `tsort`, `tclip`, `tview`, `tselect`, `tmask`, `tshift`, and grouped statistics.

### Changed

- **Breaking**: Move DimensionalData and Unitful integrations from hard dependencies to package extensions.
- **Breaking**: `tderiv` now returns plain floating-point derivatives per second for floating-point data with `Date`/`DateTime` coordinates, so Unitful is no longer required.
- **Breaking**: Remove container-target interpolation such as `tinterp(A, B::AbstractDimArray)`; pass target coordinates explicitly with `tinterp(A, times(B))`.
