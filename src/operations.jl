using DimensionalData.Lookups: Unordered

dimtype_eltype(d) = (basetypeof(d), eltype(d))
dimtype_eltype(A, ::Nothing) = dimtype_eltype(A, TimeDim)
dimtype_eltype(A, dim) = dimtype_eltype(dims(A, dim))

"""
    tsort(A; dim=nothing, rev=false)

Sort a `DimArray` `A` along the `dim` dimension.
"""
function tsort(A; dim = nothing, rev = false)
    tdim = timedim(A, dim)

    return if issorted(tdim; rev)
        DimensionalData.order(tdim) isa Unordered ?
            set(A, tdim => rev ? ReverseOrdered : ForwardOrdered) : A
    else
        time = parent(lookup(tdim))
        order = rev ? ReverseOrdered : ForwardOrdered
        sel = rebuild(tdim, sortperm(time; rev))
        set(A[sel], tdim => order)
    end
end

"""
    tclip(A, t0, t1; dim=nothing)

Clip `A` to time range `[t0, t1]` along dimension `dim`.

`dim` should be sorted before clipping (see [`tsort`](@ref)).
"""
function tclip(A, t0, t1; dim = nothing)
    Dim, T = dimtype_eltype(A, dim)
    return A[Dim(T(t0) .. T(t1))]
end

"""
    tview(A, t0, t1; dim=nothing)

View `A` in time range `[t0, t1]` along dimension `dim`.
"""
function tview(A, t0, t1; dim = nothing)
    Dim, T = dimtype_eltype(A, dim)
    return @view A[Dim(T(t0) .. T(t1))]
end

"""
    tclips(xs...; trange=common_timerange(xs...))

Clip multiple arrays to a common time range `trange`.

If `trange` is not provided, automatically finds the common time range
across all input arrays.
"""
tclips(xs::Vararg{Any, N}; trange = common_timerange(xs...)) where {N} =
    _tclips(tclip, trange, xs...)

tviews(xs::Vararg{Any, N}; trange = common_timerange(xs...)) where {N} =
    _tclips(tview, trange, xs...)

_tclips(f, ::Nothing, xs...) = throw(ArgumentError("No common time range found"))
_tclips(f, trange, xs...) =
    ntuple(i -> f(xs[i], trange...), length(xs))


"""
    tmask!(A, t0, t1; dim=nothing)
    tmask!(A, its; dim=nothing)

Mask all data values within the specified time range(s) `(t0, t1)` / `its` with NaN.
"""
function tmask!(A, t0, t1; dim = nothing)
    nan = eltype(A)(NaN)
    tview(A, t0, t1; dim) .= nan
    return A
end

function tmask!(A, its; kw...)
    for it in its
        tmask!(A, it...; kw...)
    end
    return A
end

"""
    tmask(A, args...; kwargs...)

Non-mutable version of `tmask!`. See also [`tmask!`](@ref).
"""
tmask(A, args...; kwargs...) = tmask!(copy(A), args...; kwargs...)

function tselect(times, t)
    sorted = issorted(times) ? times : sort(times)
    return sorted[searchsortednearest(sorted, t)]
end

"""
    tselect(A, t; dim=nothing)
    tselect(A, t, δt; dim=nothing)

Select the value of `A` at time nearest to `t`.

With `δt`, restricts the search to `[t-δt, t+δt]` and returns `missing` if that window is empty.
"""
function tselect(A, t, δt; dim = nothing)
    tmp = tview(A, t - δt, t + δt; dim)
    return isempty(tmp) ? missing : tselect(tmp, t; dim)
end

"""
    tshift(x, t0 = nothing; dim=nothing)

Shift the `dim` of `x` by `t0`.
"""
function tshift(x, t0 = nothing; dim = nothing)
    d = dimnum(x, dim)
    td = dims(x, d)
    times = axiskeys(x, d)
    times′ = times .- (@something t0 first(times))
    return set(x, td => times′)
end
