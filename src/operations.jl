dimtype_eltype(d) = (basetypeof(d), eltype(d))
dimtype_eltype(A, query) = dimtype_eltype(dims(A, query))
dimtype_eltype(A, ::Nothing) = dimtype_eltype(A, TimeDim)

"""
    tsort(A; query=nothing, rev=false)

Sort a `DimArray` `A` along the `query` dimension.
"""
function tsort(A; query=nothing, rev=false)
    tdim = timedim(A, query)
    issorted(tdim) ? A : begin
        time = parent(lookup(tdim))
        order = rev ? ReverseOrdered : ForwardOrdered
        sel = rebuild(tdim, sortperm(time; rev))
        set(A[sel], tdim => order)
    end
end

"""
    tclip(A, t0, t1; query=nothing, sort=false)

Clip a `Dimension` or `DimArray` `A` to a time range `[t0, t1]`.

For unordered dimensions, the dimension should be sorted before clipping (see [`tsort`](@ref)).
"""
function tclip(A::AbstractDimArray, t0, t1; query=nothing)
    Dim, T = dimtype_eltype(A, query)
    return A[Dim(T(t0) .. T(t1))]
end

"""
    tview(d, t0, t1)

View a dimension or `DimArray` in time range `[t0, t1]`.
"""
tview(d, t0, t1) = @view d[DateTime(t0)..DateTime(t1)]
function tview(da::AbstractDimArray, t0, t1; query=nothing)
    Dim, T = dimtype_eltype(da, query)
    return @view da[Dim(T(t0) .. T(t1))]
end

# Type-stable
function _tmask!(da, t0, t1, Dim, T)
    nan = NaN * oneunit(eltype(da))
    da[Dim(T(t0) .. T(t1))] .= nan
    return da
end

"""
    tmask!(da, t0, t1)
    tmask!(da, it::Interval)
    tmask!(da, its)

Mask all data values within the specified time range(s) `(t0, t1)` / `it` / `its` with NaN.
"""
tmask!(da, t0, t1; query=TimeDim) = _tmask!(da, t0, t1, dimtype_eltype(da, query)...)
function tmask!(da, its::AbstractArray; kw...)
    for it in its
        tmask!(da, it; kw...)
    end
    return da
end

function tselect!(da, t; query=nothing)
    Dim, T = dimtype_eltype(da, query)
    return da[Dim(Near(T(t)))]
end

"""
    tselect(A, t, [δt]; query=nothing)

Select the value of `A` closest to time `t` within the time range `[t-δt, t+δt]`.

Similar to `DimensionalData.Dimensions.Lookups.At` but choose the closest value and return `missing` if the time range is empty.
"""
function tselect(A, t, δt; query=nothing)
    Dim, T = dimtype_eltype(A, query)
    tmp = @views A[Dim(T(t - δt) .. T(t + δt))]
    length(tmp) == 0 ? missing : tmp[Dim(Near(T(t)))]
end

"""
    tmask(da, args...; kwargs...)

Non-mutable version of `tmask!`. See also [`tmask!`](@ref).
"""
tmask(da, args...; kwargs...) = tmask!(deepcopy(da), args...; kwargs...)

"""
    tshift(x; dim=TimeDim, t0=nothing, new_dim=nothing)

Shift the `dim` of `x` by `t0`.
"""
function tshift(x, t0=nothing; query=nothing, dim=nothing, new_dim=nothing)
    dim = @something dim dimnum(x, something(query, TimeDim))
    td = dims(x, dim)
    times = parent(lookup(td))
    t0 = @something t0 first(times)
    new_dim = @something new_dim Dim{Symbol("Time after ", t0)}
    set(x, dim => new_dim(times .- t0))
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
tclips(xs::Vararg{Any,N}; trange=common_timerange(xs...)) where {N} =
    ntuple(i -> tclip(xs[i], trange...), N)

tviews(xs::Vararg{Any,N}; trange=common_timerange(xs...)) where {N} =
    ntuple(i -> tview(xs[i], trange...), N)