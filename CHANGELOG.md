# Changelog

## Unreleased

### Changed

- **Breaking**: `cadence`, `resolution`, `samplingrate`, and `find_continuous_timeranges` take time coordinates only now; pass `times(x)` for a container.
- **Breaking**: Rename `timerange(x)` to `trange(x)`; it accepts containers only. `common_timerange(xs...)` becomes `trange(xs...)`.
- **Breaking**: `times(x)` on a plain array errors instead of returning `x`; so do `tminimum`, `tmaximum`, `targmin`, and `targmax`.
- **Breaking**: Plain arrays with more than one dimension have no default time dimension: calls without `dim`, including explicit-times forms like `tinterp(A, t, t_new)`, throw `ArgumentError`. Plain vectors are unaffected.
- **Breaking**: `mapsamples(f, A, args...)` calls `f` with the matching sample of each argument, each sliced along its own time dimension. Time is no longer passed implicitly: `mapsamples(f, A)` becomes `mapsamples(f, A, times(A))`. Keyed arguments must share the time keys of `A`.
- **Breaking**: Remove the `trange` keyword of `tclips` and `tviews`; use `tclip.(xs, Ref(tr))`.
- **Breaking**: `samplingrate(ts)` is defined without Unitful and returns hertz as a plain number.
- `times(x, dim)` is generic over `axiskeys` and `tdimnum`; container backends no longer define `times`.

### Removed

- **Breaking**: The SpaceDataModel extension; SpaceDataModel is to implement the container interface for its own variables.
- **Breaking**: `resolution`; use `cadence`.
- **Breaking**: `tview(f::Function, t0, t1)`.
- **Breaking**: `Tinterp`; external interpolators passed as `interp` receive raw time values (`Dates.value`) directly.
- **Breaking**: `time_grid`; use `tresample(A, dt)` or `range(trange(x)...; step = dt)`.
- **Breaking**: `tnorm`, `tdot`, `tcross`, `tsproj`, `tproj`, `toproj`, `tnorm_combine`; use `mapsamples(norm, B)`, `mapsamples(cross, E, B)`, `mapsamples(oproj, A, B)`, ...

## v0.2.0

### Added

- Move keyed time operations to the generic axis API, enabling `AxisKeys.KeyedArray` support for operations such as `tsort`, `tclip`, `tview`, `tselect`, `tmask`, `tshift`, and grouped statistics.

### Changed

- **Breaking**: Move DimensionalData and Unitful integrations from hard dependencies to package extensions.
- **Breaking**: `tderiv` now returns plain floating-point derivatives per second for floating-point data with `Date`/`DateTime` coordinates, so Unitful is no longer required.
- **Breaking**: Remove container-target interpolation such as `tinterp(A, B::AbstractDimArray)`; pass target coordinates explicitly with `tinterp(A, times(B))`.
