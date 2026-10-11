# Reference:
# - https://github.com/brenhinkeller/NaNStatistics.jl
# - https://github.com/JuliaSIMD/VectorizedStatistics.jl

"""
    tstat(f, x; dim=nothing)
    tstat(f, x, every; period=every, start_by=:window, dim=nothing)

Statistic `f` of each channel of `x` along time (or `dim`): over all samples, or over each window
`[start, start + period)` with starts `every` apart, keyed by window start.

`f` receives one channel's samples as a vector, empty for windows without samples.
Starts begin at the first time floored to a multiple of `every` (`start_by=:window`) or at the first time
(`start_by=:datapoint`).
"""
function tstat(f, x; dim = nothing)
    d = dimnum(x, dim)
    data = unwrap(x)
    return ndims(data) == 1 ? f(data) : _dropdim(x, map(f, eachslice(data; dims = other_dims(data, d))), d)
end

function tstat(f, x, every; period = every, start_by = :window, dim = nothing)
    d = dimnum(x, dim)
    data, t = unwrap(x), axiskeys(x, d)
    windows = _tgrid_windows(t, every, period, start_by)
    # index ranges of the windows, searched once for all channels
    ranges = collect(WindowedView{1}(axes(data, d), t, windows))
    A3 = _as3d(data, d)
    cols = [map(r -> f(view(A3, p, r, q)), ranges) for q in axes(A3, 3) for p in axes(A3, 1)]
    return rebuild(x, _stack_channels(data, d, cols), d, windows.coords)
end

tstat_doc(sym, desc = sym) = """
    $(Symbol(:t, sym))(x, [every]; period=every, start_by=:window, dim=nothing)

The $desc of `x` along time (or `dim`), ignoring NaNs: over all samples or per window, see [`tstat`](@ref).
"""

for (sym, desc) in (
        (:sum, "sum"),
        (:mean, "arithmetic mean"),
        (:median, "median"),
        (:var, "variance"),
        (:std, "standard deviation"),
        (:sem, "standard error of the mean"),
    )

    nanfunc = Symbol(:nan, sym)
    tfunc = Symbol(:t, sym)
    doc = tstat_doc(sym, desc)
    @eval @doc $doc $tfunc(x, arg...; kw...) = tstat($nanfunc, x, arg...; kw...)
end

# https://github.com/JuliaLang/julia/issues/54542
# Reduce offsets from the first element: Float64 of a raw epoch count (≈1.6e18 ns) cannot resolve sub-μs.
function tmean(x::AbstractArray{<:Dates.TimeType})
    t0 = first(x)
    return _shift(t0, mean(t -> Dates.value(t - t0), x))
end

function tmedian(x::AbstractArray{<:Dates.TimeType})
    t0 = first(x)
    return _shift(t0, median!([Dates.value(t - t0) for t in x]))
end

_shift(t0, Δ) = t0 + typeof(t0 - t0)(round(Int, Δ))
