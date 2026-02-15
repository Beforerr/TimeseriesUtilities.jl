# ============================================================================
# DimensionalData-specific overloads
# All DD-specific methods are consolidated here.
# ============================================================================

import DimensionalData as DD
using DimensionalData.Lookups: Unordered

# --- API delegation ---

dims(args...; kwargs...) = DD.dims(args...; kwargs...)
set(args...; kwargs...) = DD.set(args...; kwargs...)

function axiskeys(x::AbstractDimArray, dim)
    return unwrap(DD.dims(x, dim))
end

# --- rebuild protocol for DimensionalData ---

rebuild(x::AbstractDimArray, data) = DD.rebuild(x; data)
function rebuild(x::AbstractDimArray, data, newdims)
    return DD.rebuild(x, data, DimensionalData.format(newdims, data))
end

# --- unwrap ---

unwrap(x::AbstractDimArray) = parent(x)
unwrap(x::Dimension) = parent(lookup(x))

# --- dimnum ---

dimnum(x, query) = DimensionalData.dimnum(x, @something(query, TimeDim))

# --- times / timedim ---

"""
Returns the time indices of `x`.
"""
times(x::AbstractDimArray, args...) = lookup(timedim(x, args...))
times(x::AbstractDimStack, args...) = lookup(timedim(x, args...))

function timedim(x, query = nothing)
    query = something(query, TimeDim)
    qdim = DD.dims(x, query)
    return isnothing(qdim) ? DD.dims(x, 1) : qdim
end

# --- timerange ---

timerange(times::DimensionalData.Sampled) = timerange(parent(times))
timerange(times::Dimension) = timerange(parent(times))
timerange(x::AbstractDimArray) = timerange(times(x))

# --- reduce ---

tminimum(x::AbstractDimArray; query = nothing) = tminimum(times(x, query))
tmaximum(x::AbstractDimArray; query = nothing) = tmaximum(times(x, query))
targmin(x) = times(x)[argmin(x)]
targmax(x) = times(x)[argmax(x)]

# --- resolution ---

resolution(da::AbstractDimArray; kwargs...) = resolution(times(da); kwargs...)

# --- groupby_dynamic ---

function groupby_dynamic(x::Dimension, args...; kwargs...)
    return groupby_dynamic(parent(lookup(x)), args...; kwargs...)
end

# This is faster than `DimensionalData.format(rebuild(dim, x))` if the lookup trait keeps the same
fast_rebuild_dim(dim, x) = DD.rebuild(dim, DD.rebuild(dim.val; data = x))

# --- dimtype_eltype helper ---

dimtype_eltype(d) = (basetypeof(d), eltype(d))
dimtype_eltype(A, query) = dimtype_eltype(DD.dims(A, query))
dimtype_eltype(A, ::Nothing) = dimtype_eltype(A, TimeDim)

# ============================================================================
# Operations overloads
# ============================================================================

"""
    tsort(A; query=nothing, rev=false)

Sort a `DimArray` `A` along the `query` dimension.
"""
function tsort(A::AbstractDimArray; query = nothing, rev = false)
    tdim = timedim(A, query)

    return if issorted(tdim; rev)
        DimensionalData.order(tdim) isa Unordered ?
            set(A, tdim => rev ? ReverseOrdered : ForwardOrdered) : A
    else
        time = parent(lookup(tdim))
        order = rev ? ReverseOrdered : ForwardOrdered
        sel = DD.rebuild(tdim, sortperm(time; rev))
        set(A[sel], tdim => order)
    end
end

"""
    tclip(A::AbstractDimArray, t0, t1; query=nothing)

Clip a `DimArray` `A` to a time range `[t0, t1]`.
"""
function tclip(A::AbstractDimArray, t0, t1; query = nothing)
    Dim, T = dimtype_eltype(A, query)
    return A[Dim(T(t0) .. T(t1))]
end

"""
    tview(da::AbstractDimArray, t0, t1; query=nothing)

View a `DimArray` in time range `[t0, t1]`.
"""
function tview(da::AbstractDimArray, t0, t1; query = nothing)
    Dim, T = dimtype_eltype(da, query)
    return @view da[Dim(T(t0) .. T(t1))]
end

# Type-stable tmask! for DimArrays
function _tmask!(da, t0, t1, Dim, T)
    nan = NaN * oneunit(eltype(da))
    da[Dim(T(t0) .. T(t1))] .= nan
    return da
end

tmask!(da::AbstractDimArray, t0, t1; query = nothing) = _tmask!(da, t0, t1, dimtype_eltype(da, query)...)
function tmask!(da::AbstractDimArray, its::AbstractArray; kw...)
    for it in its
        tmask!(da, it; kw...)
    end
    return da
end

function tselect(da::AbstractDimArray, t; query = nothing)
    Dim, T = dimtype_eltype(da, query)
    return da[Dim(Near(T(t)))]
end

function tselect(A::AbstractDimArray, t, δt; query = nothing)
    Dim, T = dimtype_eltype(A, query)
    tmp = @views A[Dim(T(t - δt) .. T(t + δt))]
    return length(tmp) == 0 ? missing : tmp[Dim(Near(T(t)))]
end

"""
    tshift(x, t0=nothing; query=nothing, dim=nothing)

Shift the time dimension of a DimArray/KeyedArray by `t0`.
"""
function tshift(x::AbstractDimArray, t0 = nothing; query = nothing, dim = nothing)
    dim = @something dim dimnum(x, query)
    td = DD.dims(x, dim)
    ts = axiskeys(x, dim)
    ts′ = ts .- (@something t0 first(ts))
    return set(x, td => ts′)
end

# ============================================================================
# Stats overloads
# ============================================================================

function tstat(f, x::AbstractDimArray; dim = nothing, query = nothing)
    dim = @something dim dimnum(x, query)
    return ndims(x) == 1 ? f(x) : f(x; dim)
end

function tstat(f, x::AbstractDimArray, dt; dim = nothing, query = nothing)
    dim = @something dim dimnum(x, query)
    tdim = DD.dims(x, dim)
    out, idxs = stat1d(f, parent(x), tdim, dt, dim)
    newdims = ntuple(ndims(x)) do i
        i == dim ? fast_rebuild_dim(tdim, idxs) : DD.dims(x, i)
    end
    return rebuild(x, out, newdims)
end

function tstat(f, ds::DimStack, args...; query = nothing, dim = nothing)
    dim = @something dim dimnum(ds, query)
    return maplayers(ds) do layer
        tstat(f, layer, args...; dim)
    end
end

# ============================================================================
# Algebra overloads
# ============================================================================

function tsubtract(x::AbstractDimArray, f = nanmedian; dim = nothing, query = nothing)
    d = @something dim dimnum(x, query)
    return x .- f(parent(x); dims = d)
end

function tderiv(A::AbstractDimArray; dim = nothing, query = nothing, lazy = false)
    d = @something dim dimnum(A, query)
    tdim = DD.dims(A, d)
    f = lazy ? DiffQ : _tderiv
    out = f(parent(A), unwrap(tdim); dim = d)
    newdims = ntuple(i -> i == d ? @view(tdim[1:(end - 1)]) : DD.dims(A, i), ndims(A))
    return rebuild(A, out, newdims)
end

function tnorm_combine(x::AbstractDimArray; dim = nothing, query = nothing, name = :magnitude)
    d = @something dim dimnum(x, query)
    data = norm_combine(parent(x), d)

    odim = otherdims(x, d) |> only
    odimType = basetypeof(odim)
    new_odim = odimType(vcat(odim.val, name))
    new_dims = map(dd -> dd isa odimType ? new_odim : dd, DD.dims(x))
    return rebuild(x, data, new_dims)
end

# ============================================================================
# Interpolation overloads
# ============================================================================

function tinterp(A::AbstractDimArray, t; query = nothing, dim = nothing, kws...)
    dim = @something dim dimnum(A, query)
    out = tinterp(parent(A), unwrap(DD.dims(A, dim)), t; dim, kws...)
    return if t isa AbstractTime
        out
    else
        newdim = DD.rebuild(DD.dims(A, dim), t)
        newdims = ntuple(i -> i == dim ? newdim : DD.dims(A, i), ndims(A))
        rebuild(A, out, newdims)
    end
end

function tresample(A::AbstractDimArray, dt; query = nothing, dim = nothing, kws...)
    dim = @something dim dimnum(A, query)
    old_times = unwrap(DD.dims(A, dim))
    return tinterp(A, time_grid(old_times, dt); dim, kws...)
end

"""
    tinterp(A, B::AbstractDimArray; kws...)

Interpolate `A` to times in `B`
"""
tinterp(A, B::AbstractDimArray; kws...) = tinterp(A, times(B); kws...)

function tinterp_nans(da::AbstractDimArray; query = nothing, kwargs...)
    u = parent(da)
    dim = timedim(da, query)
    t = parent(lookup(dim))
    new_data = mapslices(u; dims = dimnum(da, dim)) do slice
        interpolate_nans(slice, t; kwargs...)
    end
    return DD.rebuild(da; data = new_data)
end

@views function tsync(A, Bs...)
    tr = common_timerange(A, Bs...)
    @assert !isnothing(tr) "No common time range found"
    A_tsync = A[Ti(Between(tr...))]
    return ntuple(1 + length(Bs)) do i
        i == 1 ? A_tsync : tinterp(Bs[i - 1], A_tsync)
    end
end

function workload_interp_setup(n = 4)
    times1 = DateTime(2020, 1, 1) + Day.(0:(n - 1))
    times2 = DateTime(2020, 1, 2) + Day.(0:(n - 1))
    times3 = DateTime(2020, 1, 1, 12) + Day.(0:(n - 2))

    da1 = DimArray(1:n, (Ti(times1),))
    da2 = DimArray(10:(10 + n - 1), (Ti(times2),))
    da3 = DimArray(hcat(5:(5 + n - 2), 8:2:(8 + 2n - 4)), (Ti(times3), Y([1, 2])))
    return da1, da2, da3
end

function workload_interp()
    da1, da2, da3 = workload_interp_setup()
    return tsync(da1, da2, da3)
end

# ============================================================================
# Groupby overloads
# ============================================================================

function tgroupby(x::AbstractDimArray, args...; dim = nothing, query = nothing, kwargs...)
    d = @something dim dimnum(x, query)
    ts = DD.dims(x, d)
    group_idx, = groupby_dynamic(ts, args...; kwargs...)
    return map(group_idx) do idx
        selectdim(x, d, idx)
    end
end

# ============================================================================
# Timeseries overloads
# ============================================================================

function smooth(da::AbstractDimArray, window; kwargs...)
    return smooth(da, Integer(div(window, resolution(da))); kwargs...)
end

function smooth(da::AbstractDimArray, window::Integer; dim = nothing, query = nothing, suffix = "_smoothed", kwargs...)
    d = @something dim dimnum(da, query)
    new_data = mapslices(parent(da); dims = d) do slice
        nanmean.(RollingWindowArrays.rolling(slice, window; kwargs...))
    end
    return DD.rebuild(da; data = new_data, name = Symbol(da.name, suffix))
end

function tfilter(da::AbstractDimArray, Wn1, Wn2 = 0.999 * samplingrate(da) / 2; designmethod = nothing)
    designmethod = something(designmethod, Butterworth(2))
    fs = samplingrate(da)
    Wn1, Wn2, fs = (Wn1, Wn2, fs) ./ 1u"Hz" .|> NoUnits
    f = digitalfilter(Bandpass(Wn1, Wn2; fs), designmethod)
    res = filtfilt(f, ustrip(parent(da)))
    return DD.rebuild(da; data = res * (da |> eltype |> unit))
end

function dropna(A::AbstractDimArray; query = nothing, dim = nothing)
    d = @something dim dimnum(A, query)
    valid_idx = vec(all(!isnan, A; dims = other_dims(A, d)))
    return selectdim(A, d, valid_idx)
end

function dropna(ds::DimStack, query = nothing)
    query = something(query, TimeDim)
    Dim, T = dimtype_eltype(ds, query)
    odims = otherdims(ds, query)

    valid_idx = mapreduce(.*, values(ds)) do A
        vec(all(!isnan, A; dims = odims))
    end

    return ds[Dim(valid_idx)]
end

# ============================================================================
# Outliers overloads
# ============================================================================

function find_outliers(A::AbstractDimArray, args...; dim = nothing, query = TimeDim, kw...)
    d = something(dim, dimnum(A, query))
    return find_outliers(parent(A), args...; dim = d, kw...)
end
