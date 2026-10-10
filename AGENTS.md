# Agent Rules

Keep generic time-series behavior in main code and package-specific details in extensions.

Prefer defining operations against common coordinate API:

- `dimnum(x, dim)`: dimension index; `dim=nothing` resolves to the time dimension.
- `axiskeys(x, dim)`: raw coordinate values for a keyed/dimensional container axis.
- `tdimnum(x)`: time dimension index, or `nothing`; backends implement this rather than `dimnum(x, nothing)`.
- `times(x, dim=nothing)`: raw time coordinate values.
- `unwrap(x)`: parent underlying array.
- `rebuild(x, data, dim, keys)`: rebuild the original container with new data and axis keys.

Container backends (extensions) implement `tdimnum`, `dimnum` for their dimension-name types,
`axiskeys`, `unwrap`, `rebuild`, and `rekey`.
