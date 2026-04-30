function resolution(times; tol = 2, check = true, f = _median)
    dts = diff(times)
    return if check
        dt, relerr = stat_relerr(dts, f)
        relerr > exp10(-tol - 1) && @warn "Time resolution is is not approximately constant (relerr ≈ $relerr)"
        dt
    else
        f(dts)
    end
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
function smooth(data, coords, window; dim = ndims(data), op = nanmean)
    length(coords) == size(data, dim) || throw(DimensionMismatch("length(coords) must match size(data, dim)"))
    issorted(coords) || throw(ArgumentError("coords must be sorted"))
    before, after = _window_offsets(window)
    return mapslices(data; dims = dim) do slice
        op.(SlidingWindow(slice, coords, before, after))
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
