"""
    tstat(f, x; dim=nothing)
    tstat(f, x, every; period=every, start_by=:window, dim=nothing)
    tstat(f, x, windows; dim=nothing)

Statistic `f` of each channel of `x` along time (or `dim`): over all samples, over each window
`[start, start + period)` with starts `every` apart, or over each of the given time ranges `windows`. Results are keyed by window start.

`f` receives one channel's samples as a vector, empty for windows without samples.
Starts begin at the first time floored to a multiple of `every` (`start_by=:window`) or at the first time
(`start_by=:datapoint`).
"""
function tstat(f, x; dim = nothing)
    return _mapslices((v, _) -> f(v), x, resolve_dim(x, dim))
end

# `f(slice, I)` of each slice of `x` along `d`, `I` indexing the other dimensions.
@inline function _mapslices(f::F, x, d) where {F}
    data = unwrap(x)
    ndims(data) == 1 && return f(data, ())
    slices = eachslice(data; dims = other_dims(data, d))
    return _dropdim(x, _map(I -> f(slices[I], Tuple(I)), CartesianIndices(slices)), d)
end

# Single entry point: avoid ambiguous methods with the `AbstractDimStack`
tstat(f, x, w; dim = nothing, kw...) = _tstat(f, x, w, resolve_dim(x, dim); kw...)

@inline function _tstat(f, x, every::Union{Number, AbstractTime}, d; period = every, start_by = :window)
    k, t, windows = _grid_windows(x, d, every, period, start_by)
    # index ranges, searched once for all channels
    ranges = collect(WindowedView{1}(axes(unwrap(x), d), t, windows))
    return _rebuild_ordered(x, _rangestat(f, unwrap(x), d, ranges), d, isnothing(k) ? nothing : windows.coords)
end

@inline function _tstat(f, x, windows, d)
    t = _checksorted(times(x, d))
    ranges = [_range_indices(t, w) for w in windows]
    return rebuild(x, _rangestat(f, unwrap(x), d, ranges), d, first.(_bounds.(windows)))
end

@inline function _rangestat(f::F, data, d, ranges) where {F}
    A3 = _as3d(data, d)
    cols = [_map(r -> f(view(A3, p, r, q)), ranges) for q in axes(A3, 3) for p in axes(A3, 1)]
    return _stack_channels(data, d, cols)
end

tstat_doc(sym, desc = sym, note = "") = """
    $(Symbol(:t, sym))(x, [every]; period=every, start_by=:window, dim=nothing)

The $desc of `x` along time (or `dim`), ignoring NaNs: over all samples or per window, see [`tstat`](@ref).
$note
"""

for (sym, f, desc, note) in (
        (:sum, :nansum, "sum", ""),
        (:mean, :_nanmean, "arithmetic mean", ""),
        (:median, :_nanmedian, "median", ""),
        (:var, :nanvar, "variance", ""),
        (:std, :nanstd, "standard deviation", ""),
        (:sem, :nansem, "standard error of the mean", ""),
        (:maximum, :nanmaximum, "maximum", "Empty windows of integer data give `typemin`."),
        (:minimum, :nanminimum, "minimum", "Empty windows of integer data give `typemax`."),
    )

    tfunc = Symbol(:t, sym)
    doc = tstat_doc(sym, desc, note)
    @eval @doc $doc $tfunc(x, arg...; kw...) = tstat($f, x, arg...; kw...)
end

_nanmean(x) = nanmean(x)
_nanmedian(x) = nanmedian(x)

# https://github.com/JuliaLang/julia/issues/54542
# Reduce offsets from the first element: Float64 of a raw epoch count (≈1.6e18 ns) cannot resolve sub-μs.
_nanmean(x::AbstractArray{<:Dates.TimeType}) = _offset_stat(mean, x)
_nanmedian(x::AbstractArray{<:Dates.TimeType}) = _offset_stat(median, x)

function _offset_stat(f, x)
    t0 = first(x)
    return t0 + typeof(t0 - t0)(round(Int, f(Dates.value(t - t0) for t in x)))
end
