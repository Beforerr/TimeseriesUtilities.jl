"""
    tsort(A; dim=nothing, rev=false)

Sort `A` along dimension `dim`.
"""
function tsort(A; dim = nothing, rev = false)
    d = resolve_dim(A, dim)
    keys = axiskeys(A, d)
    return issorted(keys; rev) ? A : begin
            sorted = A[_selectors(A, d, sortperm(keys; rev))...]
            sorted_axis(sorted, d; rev)
        end
end

"""
    tclip(A, t0, t1; dim=nothing)
    tclip(A, tr; dim=nothing)

Clip `A` to time range `[t0, t1]`, or to `tr`, along dimension `dim`.
`tr` is a `(t0, t1)` tuple or an interval such as `IntervalSets.Interval{:closed, :open}(t0, t1)`.

`dim` should be sorted before clipping (see [`tsort`](@ref)).
"""
tclip(A, t0, t1; dim = nothing) = tclip(A, (t0, t1); dim)

function tclip(A, tr; dim = nothing)
    d = resolve_dim(A, dim)
    idx = _search_range(axiskeys(A, d), tr)
    return A[_selectors(A, d, idx)...]
end

"""
    tview(A, t0, t1; dim=nothing)
    tview(A, tr; dim=nothing)

View `A` in time range `[t0, t1]`, or in `tr` (see [`tclip`](@ref)), along dimension `dim`.
"""
tview(A, t0, t1; dim = nothing) = tview(A, (t0, t1); dim)

function tview(A, tr; dim = nothing)
    d = resolve_dim(A, dim)
    idx = _search_range(axiskeys(A, d), tr)
    return @view A[_selectors(A, d, idx)...]
end

function _selectors(A, d, selector)
    return ntuple(i -> i == d ? selector : Colon(), ndims(A))
end

function _search_range(keys, tr)
    issorted(keys) || throw(ArgumentError("axis keys must be sorted before clipping"))
    return _range_indices(keys, tr)
end

_range_indices(keys, (t0, t1)) = searchsortedfirst(keys, t0):searchsortedlast(keys, t1)


"""
    tmask!(A, t0, t1; dim=nothing)
    tmask!(A, trs; dim=nothing)

Set all data values within time range `[t0, t1]`, or within each range of `trs` to NaN.
"""
tmask!(A, t0, t1; dim = nothing) = tmask!(A, ((t0, t1),); dim)

function tmask!(A, trs; dim = nothing)
    d = resolve_dim(A, dim)
    data, keys = unwrap(A), axiskeys(A, d)
    _check_holds_nan(data, :tmask)
    for tr in trs
        view(data, _selectors(data, d, _search_range(keys, tr))...) .= eltype(data)(NaN)
    end
    return A
end

"""
    tmask(A, args...; kwargs...)

Copy of `A` (promoted to float) masked by [`tmask!`](@ref).
"""
tmask(A, args...; kwargs...) = tmask!(_float_copy(A), args...; kwargs...)

function tselect(A, t; dim = nothing)
    d = resolve_dim(A, dim)
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
    return isempty(unwrap(tmp)) ? missing : tselect(tmp, t; dim)
end

"""
    tshift(x, t0 = nothing; dim=nothing)

Shift the `dim` of `x` by `t0`.
"""
function tshift(x, t0 = nothing; dim = nothing)
    d = resolve_dim(x, dim)
    times = axiskeys(x, d)
    times′ = times .- (@something t0 first(times))
    return rekey(x, d, times′)
end
