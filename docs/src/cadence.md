# Sampling cadence

`cadence` (alias `resolution`) estimates the fundamental sampling step `Δt` of
a time series. Real timestamps rarely arrive on a clean grid — they may have:

- **gaps** at multiples of `Δt` (`k·Δt`),
- **jitter** around the nominal step (`Δt ± ε`),
- **duplicates** or out-of-order samples,
- **outliers** (single huge gaps from dropouts),
- **type quirks** (`Dates.Period` vs `Float64` timestamps).

A plain `median(diff(t))` collapses on any of these once gaps dominate. The
algorithm below recovers the true `Δt` instead.

## Algorithm

1. Diff and drop non-positive.
2. IQR-trim the top tail (`q75 + 3·IQR`) to remove huge isolated gaps.
3. Cluster sorted diffs with relative tolerance (`cluster_rtol ≈ 5%`),
   using the running mean as the center.
4. Pick the **smallest strong** cluster — not the most populous. Gaps
   produce `2Δt`, `3Δt`, … peaks; the true `Δt` is the smallest dense one.
   "Strong" means count ≥ `strong_frac · max_count` and ≥ `min_count`.
5. Refine: median of diffs within `rtol` of that center.
6. Validate: count gaps (`> 1.5·Δt`) and non-integer-multiple intervals;
   warn when present (suppress with `check = false`).

## Key choices

- **Two-tier `rtol`.** Clustering tolerance (~5%) wider than precision
  tolerance — must absorb jitter but stay below 50% so `2Δt` separates
  from `Δt`.
- **Median, not mean,** for refinement — outlier-safe.
- **Smallest strong, not plurality.** Under heavy dropout, `2Δt` or `3Δt`
  can outnumber `Δt`; plurality picks the wrong fundamental.

## Example

```@example cadence
using TimeseriesUtilities
using Dates

# Regular cadence with isolated burst gap (10 000× Δt)
t = vcat(Millisecond.(1:50), Millisecond.(10_001:10_050))
cadence(t; check = false)
```

```@example cadence
# Majority-gap dropout: only 3 of every 7 samples kept.
# Diffs are dominated by 3Δt, but the fundamental is still 1Δt.
kept = Millisecond.(vcat([[7k + 1, 7k + 2, 7k + 5] for k in 0:19]...))
cadence(kept; check = false)
```

```@example cadence
# ±3% jitter on 1 s cadence
t_jitter = (1:200) .+ 0.03 .* sin.(1:200)
cadence(t_jitter; check = false)
```

```@docs
cadence
resolution
```


## Rejected alternatives

- `median(diff)` — fails when gaps dominate.
- Histogram mode — bin width depends on `Δt`.
- Numerical greatest common divisor — breaks under jitter.