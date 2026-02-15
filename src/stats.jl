# Reference:
# - [NaNStatistics.jl](https://github.com/brenhinkeller/NaNStatistics.jl)
# - [VectorizedStatistics.jl](https://github.com/JuliaSIMD/VectorizedStatistics.jl)
# - [Average of Dates · Issue · JuliaLang/julia](https://github.com/JuliaLang/julia/issues/54542)

@inline function stat1d(f, x, index, dt, dim)
    group_idx, idxs = groupby_dynamic(index, dt)
    out = mapslices(x; dims = dim) do slice
        map(group_idx) do idx
            f(view(slice, idx))
        end
    end
    return out, idxs
end


"""
    tstat(f, x, [dt]; dim=1)

Calculate the statistic `f` of `x` along the `dim` dimension, optionally grouped by `dt`.

See also: [`groupby_dynamic`](@ref)
"""
function tstat end

function tstat(f, x; dim = 1)
    return ndims(x) == 1 ? f(x) : f(x; dim)
end

function tstat(f, x, index, dt; dim = 1)
    out, idxs = stat1d(f, x, index, dt, dim)
    return out, idxs
end


tstat_doc(sym, desc = sym) = """
    $(Symbol(:t, sym))(x, [dt]; dim=nothing, query=nothing)

Calculate the $desc of `x` along the `dim` dimension, optionally grouped by `dt`.

It returns a value if `x` is a vector along the `dim` dimension, otherwise returns an array with the specified dimension reduced.

If `dim` is not specified, it defaults to `1` (or the `query` dimension when using DimensionalData).
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

# https://github.com/JuliaLang/julia/issues/54542"
tmean(vec::AbstractArray{DateTime}) = convert(DateTime, Millisecond(round(nanmean(Dates.value.(vec)))))
