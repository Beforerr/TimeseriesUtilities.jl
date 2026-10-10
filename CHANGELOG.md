# Changelog

## Unreleased

### Changed

- **Breaking**: `cadence`, `resolution`, `samplingrate`, and `find_continuous_timeranges` take time coordinates only now; pass `times(x)` for a container.
- **Breaking**: Rename `timerange(x)` to `trange(x)`; it accepts containers only. `common_timerange(xs...)` becomes `trange(xs...)`.
- **Breaking**: `times(x)` on a plain array errors instead of returning `x`; so do `tminimum`, `tmaximum`, `targmin`, and `targmax`.
- **Breaking**: Plain arrays with more than one dimension have no default time dimension: calls without `dim`, including explicit-times forms like `tinterp(A, t, t_new)`, throw `ArgumentError`. Plain vectors are unaffected.
- `times(x, dim)` is generic over `axiskeys` and `tdimnum`; container backends no longer define `times`.
- **Breaking**: Remove the `trange` keyword of `tclips` and `tviews`; use `tclip.(xs, Ref(tr))`.

### Removed

- **Breaking**: `time_grid`; use `tresample(A, dt)` or `range(trange(x)...; step = dt)`.

## v0.2.0

### Added

- Move keyed time operations to the generic axis API, enabling `AxisKeys.KeyedArray` support for operations such as `tsort`, `tclip`, `tview`, `tselect`, `tmask`, `tshift`, and grouped statistics.

### Changed

- **Breaking**: Move DimensionalData and Unitful integrations from hard dependencies to package extensions.
- **Breaking**: `tderiv` now returns plain floating-point derivatives per second for floating-point data with `Date`/`DateTime` coordinates, so Unitful is no longer required.
- **Breaking**: Remove container-target interpolation such as `tinterp(A, B::AbstractDimArray)`; pass target coordinates explicitly with `tinterp(A, times(B))`.
