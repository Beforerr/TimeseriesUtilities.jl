# Generic operations on (data, times) pairs
# DimensionalData-specific overloads are in DimensionalData.jl

# --- Generic helpers ---

"""Find indices where t0 <= times <= t1"""
function _time_indices(ts, t0, t1)
    if issorted(ts)
        i0 = searchsortedfirst(ts, t0)
        i1 = searchsortedlast(ts, t1)
        return i0:i1
    else
        return findall(t -> t0 <= t <= t1, ts)
    end
end

# --- tsort ---

"""
    tsort(data, ts; dim=1, rev=false)
    tsort(A; query=nothing, rev=false)

Sort data along the time dimension. The generic form takes a data array and time vector.
"""
function tsort(data::AbstractArray, ts::AbstractVector; dim = 1, rev = false)
    perm = sortperm(ts; rev)
    return selectdim(data, dim, perm), ts[perm]
end

# --- tclip ---

"""
    tclip(data, ts, t0, t1; dim=1)
    tclip(A, t0, t1; query=nothing)

Clip data to a time range `[t0, t1]`. The generic form takes a data array and time vector.

For unordered dimensions, the dimension should be sorted before clipping (see [`tsort`](@ref)).
"""
function tclip(data::AbstractArray, ts::AbstractVector, t0, t1; dim = 1)
    idx = _time_indices(ts, t0, t1)
    return selectdim(data, dim, idx), ts[idx]
end

# --- tview ---

"""
    tview(data, ts, t0, t1; dim=1)
    tview(d, t0, t1)

View data in time range `[t0, t1]`. The generic form takes a data array and time vector.
"""
function tview(data::AbstractArray, ts::AbstractVector, t0, t1; dim = 1)
    idx = _time_indices(ts, t0, t1)
    return Base.selectdim(Base.view(data, ntuple(Returns(Colon()), ndims(data))...), dim, idx), view(ts, idx)
end

# Fallback for plain vector of times
tview(d, t0, t1) = @view d[DateTime(t0) .. DateTime(t1)]

# --- tmask! ---

"""
    tmask!(data, ts, t0, t1; dim=1)
    tmask!(da, t0, t1)
    tmask!(da, it::Interval)
    tmask!(da, its)

Mask all data values within the specified time range(s) `(t0, t1)` / `it` / `its` with NaN.
"""
function tmask!(data::AbstractArray, ts::AbstractVector, t0, t1; dim = 1)
    idx = _time_indices(ts, t0, t1)
    nan = NaN * oneunit(eltype(data))
    selectdim(data, dim, idx) .= nan
    return data
end

function tmask!(data::AbstractArray, ts::AbstractVector, its::AbstractArray; kw...)
    for it in its
        tmask!(data, ts, it...; kw...)
    end
    return data
end

"""
    tmask(data, ts, args...; kwargs...)

Non-mutable version of `tmask!`. See also [`tmask!`](@ref).
"""
tmask(da, args...; kwargs...) = tmask!(deepcopy(da), args...; kwargs...)

# --- tselect ---

"""
    tselect(times, t)
    tselect(data, ts, t; dim=1)
    tselect(data, ts, t, δt; dim=1)

Select the value of data closest to time `t`, optionally within the time range `[t-δt, t+δt]`.
Returns `missing` if the time range is empty (when `δt` is specified).
"""
function tselect(times, t)
    idx = issorted(times) ? searchsortednearest(times, t) :
        searchsortednearest(sort(times), t)
    return times[idx]
end

function tselect(data::AbstractArray, ts::AbstractVector, t; dim = 1)
    idx = searchsortednearest(issorted(ts) ? ts : sort(ts), t)
    return selectdim(data, dim, idx)
end

function tselect(data::AbstractArray, ts::AbstractVector, t, δt; dim = 1)
    idx = _time_indices(ts, t - δt, t + δt)
    if isempty(idx)
        return missing
    else
        sub_ts = ts[idx]
        nearest = searchsortednearest(sub_ts, t)
        return selectdim(data, dim, idx[nearest])
    end
end

# --- tshift ---

"""
    tshift(data, ts, t0=nothing; dim=1)
    tshift(x, t0=nothing; dim=nothing)

Shift the time dimension of data by `t0`. If `t0` is not specified, shifts to start at zero.
The generic form returns `(data, new_times)`.
"""
function tshift(data::AbstractArray, ts::AbstractVector, t0 = nothing; dim = 1)
    new_ts = ts .- (@something t0 first(ts))
    return data, new_ts
end

for f in (:tclip, :tview, :tmask!, :tmask)
    @eval $f(da, trange; kwargs...) = $f(da, trange...; kwargs...)
end

for f in (:tclip, :tview, :tmask!, :tmask, :tshift)
    @eval $f(args...; kwargs...) = da -> $f(da, args...; kwargs...)
end


"""
    tclips(xs...; trange=common_timerange(xs...))

Clip multiple arrays to a common time range `trange`.

If `trange` is not provided, automatically finds the common time range
across all input arrays.
"""
tclips(xs::Vararg{Any, N}; trange = common_timerange(xs...)) where {N} =
    ntuple(i -> tclip(xs[i], trange...), N)

tviews(xs::Vararg{Any, N}; trange = common_timerange(xs...)) where {N} =
    ntuple(i -> tview(xs[i], trange...), N)
