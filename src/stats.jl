# Reference:
# - https://github.com/brenhinkeller/NaNStatistics.jl
# - https://github.com/JuliaSIMD/VectorizedStatistics.jl

@inline function stat1d(f, x, index, dt, dim)
    group_idx, coords = groupby_dynamic(index, dt)
    out = mapslices(x; dims = dim) do slice
        map(group_idx) do idx
            f(view(slice, idx))
        end
    end
    return out, coords
end


"""
    tstat(f, x, [dt]; dim = nothing)

Calculate the statistic `f` of `x` along the `dim` dimension, optionally grouped by `dt`.

See also: [`groupby_dynamic`](@ref)
"""
function tstat end

function tstat(f, x; dim = nothing)
    d = dimnum(x, dim)
    return ndims(x) == 1 ? f(x) : f(x; dim = d)
end

function tstat(f, x, dt; dim = nothing)
    d = dimnum(x, dim)
    out, s = stat1d(f, unwrap(x), axiskeys(x, d), dt, d)
    return rebuild(x, out, d, s)
end

tstat_doc(sym, desc = sym) = """
    $(Symbol(:t, sym))(x, [dt]; dim=nothing)

Calculate the $desc of `x` along dimension `dim`, optionally grouped by `dt`.

`dim` accepts integer index, dimension type/instance, or `nothing` (defaults to the time dimension).
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
