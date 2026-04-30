"""
    tsort(A; dim=nothing, rev=false)

Sort `A` along dimension `dim`.
"""
function tsort(A; dim = nothing, rev = false)
    d = dimnum(A, dim)
    keys = axiskeys(A, d)
    return issorted(keys; rev) ? A : begin
            sorted = A[_selectors(A, d, sortperm(keys; rev))...]
            sorted_axis(sorted, d; rev)
        end
end

"""
    tclip(A, t0, t1; dim=nothing)

Clip `A` to time range `[t0, t1]` along dimension `dim`.

`dim` should be sorted before clipping (see [`tsort`](@ref)).
"""
function tclip(A, t0, t1; dim = nothing)
    d = dimnum(A, dim)
    idx = _search_range(axiskeys(A, d), t0, t1)
    return A[_selectors(A, d, idx)...]
end

"""
    tview(A, t0, t1; dim=nothing)

View `A` in time range `[t0, t1]` along dimension `dim`.
"""
function tview(A, t0, t1; dim = nothing)
    d = dimnum(A, dim)
    idx = _search_range(axiskeys(A, d), t0, t1)
    return @view A[_selectors(A, d, idx)...]
end

function _selectors(A, d, selector)
    return ntuple(i -> i == d ? selector : Colon(), ndims(A))
end

function _search_range(keys, t0, t1)
    issorted(keys) || throw(ArgumentError("axis keys must be sorted before clipping"))
    i0 = searchsortedfirst(keys, t0)
    i1 = searchsortedlast(keys, t1)
    return i0:i1
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

function tselect(A, t; dim = nothing)
    d = dimnum(A, dim)
    keys = axiskeys(A, d)
    i = if issorted(keys)
        searchsortednearest(keys, t)
    else
        order = sortperm(keys)
        order[searchsortednearest(view(keys, order), t)]
    end
    return A[_selectors(A, d, i)...]
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
