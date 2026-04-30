# Changelog

## v0.2.0

### Added

- Move keyed time operations to the generic axis API, enabling `AxisKeys.KeyedArray` support for operations such as `tsort`, `tclip`, `tview`, `tselect`, `tmask`, `tshift`, and grouped statistics.

### Changed

- **Breaking**: Move DimensionalData and Unitful integrations from hard dependencies to package extensions.
- **Breaking**: `tderiv` now returns plain floating-point derivatives per second for floating-point data with `Date`/`DateTime` coordinates, so Unitful is no longer required.
- **Breaking**: Remove container-target interpolation such as `tinterp(A, B::AbstractDimArray)`; pass target coordinates explicitly with `tinterp(A, times(B))`.
