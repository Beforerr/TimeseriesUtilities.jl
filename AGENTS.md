# Agent Rules

Keep generic time-series behavior in main code and package-specific details in extensions.

Prefer defining operations against common coordinate API:

- `axiskeys(x, dim)`: raw coordinate values for a keyed/dimensional container axis.
- `dims(x, dim)`: dimension identity/metadata, not coordinate values.
- `times(x; dim=nothing)` or `times(x, dim)`: raw time coordinate values.
- `unwrap(x)`: parent underlying array.
- `rebuild_axis(x, data, dim, keys)`: rebuild the original container with new data and axis keys.
