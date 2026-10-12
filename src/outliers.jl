# - [Data Preprocessing - MATLAB](https://www.mathworks.com/help/matlab/preprocessing-data.html?s_tid=CRUX_lftnav)
# - https://www.mathworks.com/help/signal/smoothing-and-denoising.html
# - https://www.mathworks.com/help/matlab/ref/cleanoutlierdata.html

# References
# - [HampelOutliers.jl](https://github.com/tobydriscoll/HampelOutliers.jl)
# - https://www.mathworks.com/help/signal/ref/hampel.html
# - https://github.com/melizalab/quickspikes
# - https://halleysfifthinc.github.io/Peaks.jl/stable/
# - https://pyspedas.readthedocs.io/en/latest/_modules/pytplot/tplot_math/clean_spikes.html
# - https://stackoverflow.com/questions/37556487/remove-spikes-from-signal-in-python
# - https://medium.com/towards-data-science/removing-spikes-from-raman-spectra-a-step-by-step-guide-with-python-b6fd90e8ea77

"""
    outlier_mask(A, [window]; threshold=3, dim=nothing)

Boolean mask (like `A`) of the values of each channel of `A` more than `threshold` scaled median absolute deviations
(MAD) from the median along time (or `dim`); NaNs are never outliers.
The statistics are global, or local over a moving `window` around each sample as in [`smooth`](@ref)
(counted in samples for a plain array).

See also: [`despike`](@ref), [isoutlier - MATLAB](https://www.mathworks.com/help/matlab/ref/isoutlier.html)
"""
function outlier_mask(A, window = nothing; threshold = 3, dim = nothing)
    d = dimnum(A, dim)
    data = unwrap(A)
    mask = similar(data, Bool)
    windows = isnothing(window) ? nothing : _point_windows(A, d, window)
    _outlier_mask!(_as3d(mask, d), _as3d(data, d), windows, threshold)
    return rebuild(A, mask)
end

function _outlier_mask!(M3, A3, windows, k)
    sorted = Vector{float(eltype(A3))}()
    _foreach_pq(A3) do p, q
        x, m = view(A3, p, :, q), view(M3, p, :, q)
        if isnothing(windows)
            sort!(filter!(!isnan, copyto!(resize!(sorted, length(x)), x)))
            med, σ = _median_mad(sorted)
            m .= abs.(x .- med) .> k * σ
        else
            _rolling_outlier_mask!(m, x, windows, k, empty!(sorted))
        end
    end
    return M3
end

# Window ends only advance along sorted coordinates, so `lo` and `hi` move forward and `sorted` (the non-NaN
# `x[lo:hi]`) is updated by inserting the samples that entered and deleting those that left, whatever the spacing.
function _rolling_outlier_mask!(m, x, windows, k, sorted)
    c = windows.coords
    lo, hi = 1, 0
    for (i, (t_start, t_stop)) in enumerate(windows)
        while hi < length(c) && c[hi + 1] < t_stop
            hi += 1
            isnan(x[hi]) || insert!(sorted, searchsortedfirst(sorted, x[hi]), x[hi])
        end
        while lo <= hi && c[lo] < t_start
            isnan(x[lo]) || deleteat!(sorted, searchsortedfirst(sorted, x[lo]))
            lo += 1
        end
        med, σ = _median_mad(sorted)
        m[i] = abs(x[i] - med) > k * σ
    end
    return m
end

# Median and scaled MAD of sorted `s`.
function _median_mad(s)
    n = length(s)
    n == 0 && return eltype(s)(NaN), eltype(s)(NaN)
    med = _middle(i -> s[i], n)
    return med, 1.4826 * _middle(k -> _kth_deviation(s, med, k), n)
end

_middle(f, n) = (f((n + 1) ÷ 2) + f(n ÷ 2 + 1)) / 2

# k-th smallest `|s[i] - med|`: the k values of sorted `s` nearest `med` are a block `s[a:a+k-1]`;
# `a` is the first start whose left end is no farther than the value after the block.
function _kth_deviation(s, med, k)
    lo, hi = 1, length(s) - k + 1
    while lo < hi
        a = (lo + hi) ÷ 2
        med - s[a] > s[a + k] - med ? (lo = a + 1) : (hi = a)
    end
    return max(med - s[lo], s[lo + k - 1] - med)
end

"""
    despike(A, [window]; threshold=3, dim=nothing)
    despike!(A, [window]; threshold=3, dim=nothing)

Set the [`outlier_mask`](@ref) values of `A` to NaN.
"""
function despike(A, args...; kw...)
    data = unwrap(A)
    return despike!(rebuild(A, copyto!(similar(data, float(eltype(data))), data)), args...; kw...)
end

function despike!(A, args...; kw...)
    data = unwrap(A)
    eltype(data) <: AbstractFloat ||
        throw(ArgumentError("despike! sets outliers to NaN, which an array of $(eltype(data)) cannot hold; use despike"))
    data[unwrap(outlier_mask(A, args...; kw...))] .= NaN
    return A
end
