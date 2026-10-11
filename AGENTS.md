# Agent Rules

Keep generic time-series behavior in main code and package-specific details in extensions.

- See `src/api.jl` for backend hooks and their contract; `Timed` (`src/timed.jl`) is the reference backend for plain arrays.
- Write operations as kernels on `(data, t, d)` (plain data, time keys, dimension index) and lift with `_ontime`, or unwrap → kernel → `rebuild` when the output keys need more.
- Every generic operation belongs in the `ops` list of `test/backends.jl`, which checks all backends against `Timed`.
