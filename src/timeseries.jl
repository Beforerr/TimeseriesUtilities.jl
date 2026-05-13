"""
    cadence(x; tol=2, rtol=exp10(-tol-1), check=true)

Robust sampling cadence of `times(x)` via IQR-trimmed modal clustering on
positive diffs, refined by the median of values near the modal center.
With `check=true`, warns on gaps or non-integer-multiple intervals.
"""
cadence(x::AbstractRange; kwargs...) = step(x)

function cadence(x; tol = 2, rtol = exp10(-tol - 1), check = true)
    ts = times(x)
    length(ts) > 1 || throw(ArgumentError("at least two time points are required"))
    T = typeof(ts[2] - ts[1])

    pos = _positive_diffs(ts)
    isempty(pos) && throw(ArgumentError("no positive time differences found"))
    sort!(pos)

    cval = _modal_cadence(pos, rtol)
    check && _check_regular(pos, cval, rtol)
    return T <: AbstractFloat ? cval : T(round(Int, cval))
end

const resolution = cadence

function _positive_diffs(ts)
    @inline _float(x) = Float64(x)
    @inline _float(x::Union{Dates.Period, Dates.AbstractTime}) = Float64(Dates.value(x))

    n = length(ts) - 1
    out = Vector{Float64}(undef, n)
    any_bad = false
    @inbounds @simd for i in 1:n
        v = _float(ts[i + 1] - ts[i])
        out[i] = v
        any_bad |= (v <= 0.0)
    end
    any_bad && filter!(>(0.0), out)
    return out
end

# Pre-sorted positive diffs → fundamental cadence.
# Clustering rtol independent of precision rtol: wide enough to absorb jitter,
# narrow enough that 2Δt stays a separate cluster from Δt.
function _modal_cadence(
        pos, rtol;
        cluster_rtol = 0.05, strong_frac = 0.25, min_count = 3
    )
    n = length(pos)
    q25 = _quantile_sorted(pos, 0.25)
    q75 = _quantile_sorted(pos, 0.75)
    j = clamp(searchsortedlast(pos, q75 + 3 * (q75 - q25)), 1, n)

    p1, pj = pos[1], pos[j]
    # single-cluster check: all values within ±cluster_rtol of midpoint.
    best = if pj - p1 <= cluster_rtol * (p1 + pj)
        _median_sorted(pos, 1, j)
    else
        _smallest_strong_center(pos, j, cluster_rtol, strong_frac, min_count)
    end

    # Refine via median of sorted slice within ±rtol of best.
    i_lo = searchsortedfirst(pos, best * (1 - rtol))
    i_hi = searchsortedlast(pos, best * (1 + rtol))
    return i_lo <= i_hi ? _median_sorted(pos, i_lo, i_hi) : best
end

# Walk sorted pos[1:j], form running-mean clusters, return smallest center
# whose count clears strong threshold. Centers come out sorted, so findfirst suffices.
function _smallest_strong_center(pos, j, cluster_rtol, strong_frac, min_count)
    centers = Float64[]
    counts = Int[]
    maxk = 0
    @inbounds begin
        c = pos[1]; k = 1
        for i in 2:j
            v = pos[i]
            if v - c <= cluster_rtol * c
                k += 1
                c += (v - c) / k
            else
                push!(centers, c); push!(counts, k)
                k > maxk && (maxk = k)
                c = v; k = 1
            end
        end
        push!(centers, c); push!(counts, k)
        k > maxk && (maxk = k)
    end
    thresh = max(min_count, ceil(Int, strong_frac * maxk))
    idx = findfirst(>=(thresh), counts)
    return isnothing(idx) ? centers[1] : centers[idx]
end

function _check_regular(pos, cval, rtol)
    gap_thresh = 1.5 * cval
    bad_thresh = rtol * cval
    n_gaps = 0
    n_bad = 0
    @inbounds for v in pos
        v > gap_thresh && (n_gaps += 1)
        abs(v - round(v / cval) * cval) > bad_thresh && (n_bad += 1)
    end
    if n_gaps > 0 || n_bad > 0
        parts = String[]
        n_gaps > 0 && push!(parts, "$n_gaps gap(s)")
        n_bad > 0 && push!(parts, "$n_bad non-integer-multiple interval(s)")
        @warn "Time resolution is not approximately constant ($(join(parts, ", ")))"
    end
    return nothing
end

"""
    smooth(data, times, window; dim=ndims(data))
    smooth(data, window; dim=ndims(data))

Smooths a time series by computing a moving average over a sliding window.
Edge windows are truncated, so the output has the same size as the input.

Each window covers the half-open interval `[coord - before, coord + after)`.
A scalar `window` is interpreted as a coordinate span along the smoothed axis.

# Arguments
- `dim=ndims(data)`: Dimension along which to perform smoothing
- `op=nanmean`: Function used to aggregate each window
"""
@inline function smooth(data, coords, window; dim = ndims(data), op = nanmean)
    length(coords) == size(data, dim) || throw(DimensionMismatch("length(coords) must match size(data, dim)"))
    issorted(coords) || throw(ArgumentError("coords must be sorted"))
    before, after = _window_offsets(window)
    windows = PointWindows(coords, before, after)
    return mapslices(data; dims = dim) do slice
        op.(WindowedView{1}(slice, coords, windows))
    end
end

function smooth(data, window; dim = ndims(data), kw...)
    return smooth(data, axes(data, dim), window; dim, kw...)
end


function _window_offsets(window::Tuple)
    @assert length(window) == 2
    return window[1], window[2]
end

_window_offsets(window) = _half(window), _half(window)
_half(window) = window / 2
_half(window::Period) = Millisecond(window) / 2

"""
    dropna(A; dim=nothing)

Remove slices containing NaN values along dimension `dim`.
"""
@inline function dropna(A; dim = nothing)
    d = dimnum(A, dim)
    idxs = vec(all(!isnan, A; dims = other_dims(A, d)))
    return selectdim(A, d, idxs)
end
