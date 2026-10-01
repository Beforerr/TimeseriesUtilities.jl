# https://github.com/brendanjohnharris/TimeseriesTools.jl/blob/main/ext/DataInterpolationsExt.jl
# https://github.com/JuliaMath/Interpolations.jl
# https://discourse.julialang.org/t/interpolating-along-a-single-dimension-of-a-multi-dimensional-array-for-particular-points/29308/3

struct LinearInterpolation{U, T, E}
    u::U
    t::T
    extrapolation::E
end

function check_interp_args(u, t)
    issorted(t) || throw(ArgumentError("t must be sorted"))
    length(u) == length(t) || throw(DimensionMismatch("u and t must have the same length"))
    return !isempty(t) || throw(ArgumentError("at least one interpolation point is required"))
end

function LinearInterpolation(u, t; extrapolation = false, check = true)
    check && check_interp_args(u, t)
    return LinearInterpolation(u, t, extrapolation)
end

function _interp_segment(t, x, extrapolation)
    n = length(t)
    if n == 1
        (x == only(t) || extrapolation) && return 1
        throw(DomainError(x, "interpolation point is outside the time range"))
    end

    i = searchsortedlast(t, x)
    return if i == 0
        extrapolation ? 1 : throw(DomainError(x, "interpolation point is before the time range"))
    elseif i >= n
        if x == t[end]
            return n - 1
        end
        extrapolation ? n - 1 : throw(DomainError(x, "interpolation point is after the time range"))
    else
        i
    end
end

# Sorted queries land at or after the previous segment: gallop forward from `hint`
# (O(log gap) instead of O(log n)), then bisect the bracket.
@inline function _interp_segment(t, x, extrapolation, hint)
    n = length(t)
    @inbounds if 1 <= hint < n && t[hint] <= x < t[n]
        x < t[hint + 1] && return hint
        lo, step = hint + 1, 1
        hi = lo + 1
        while t[hi] <= x
            lo, step = hi, 2step
            hi = min(lo + step, n)
        end
        while hi - lo > 1
            m = (lo + hi) >>> 1
            t[m] <= x ? (lo = m) : (hi = m)
        end
        return lo
    end
    return _interp_segment(t, x, extrapolation)
end

function (interp::LinearInterpolation)(x)
    if length(interp.t) == 1
        (x == only(interp.t) || interp.extrapolation) && return only(interp.u)
        throw(DomainError(x, "interpolation point is outside the time range"))
    end
    i = _interp_segment(interp.t, x, interp.extrapolation)
    t0 = interp.t[i]
    t1 = interp.t[i + 1]
    u0 = interp.u[i]
    u1 = interp.u[i + 1]
    return u0 + (x - t0) / (t1 - t0) * (u1 - u0)
end

"""
    tinterp(A, old_times, new_times; interp=LinearInterpolation)

Interpolate time series `A` at new time points `new_times`.

The `interp` constructor must accept `(u, t; kws...)` and return a callable object.
Its syntax is compatible with `DataInterpolations.jl`.

# Examples

```julia
# Interpolate at a single time point
tinterp(time_series, DateTime("2023-01-01T12:00:00"))

# Interpolate at multiple time points using cubic spline interpolation
new_times = DateTime("2023-01-01"):Hour(1):DateTime("2023-01-02")
tinterp(time_series, new_times; interp = CubicSpline)
```
"""
function tinterp(A, old_times, new_times; interp = LinearInterpolation, dim = ndims(A), kws...)
    new_times isa AbstractArray && return _tinterp_nd(interp, A, old_times, new_times, dim; kws...)
    out = _tinterp_nd(interp, A, old_times, [new_times], dim; kws...)
    return ndims(A) == 1 ? only(out) : out
end

# Fast path for built-in LinearInterpolation: one pass over all channels per query point.
@inline _tinterp_nd(::Type{<:LinearInterpolation}, A, old_times, new_times, d; extrapolation = false, kws...) =
    _tinterp_linear_nd(A, old_times, new_times, d, extrapolation)

function _tinterp_linear_nd(A, old_times, new_times, d, extrapolation)
    out_sz = ntuple(i -> i == d ? length(new_times) : size(A, i), ndims(A))
    out = similar(A, float(eltype(A)), out_sz)
    _tinterp_linear3!(_as3d(out, d), _as3d(A, d), old_times, new_times, extrapolation)
    return out
end

function _tinterp_linear3!(out3, A3, old_times, new_times, extrapolation)
    hint = 0
    @inbounds for (j, x) in enumerate(new_times)
        i = hint = _interp_segment(old_times, x, extrapolation, hint)
        α = (x - old_times[i]) / (old_times[i + 1] - old_times[i])
        β = 1 - α
        _foreach_pq(A3) do p, q
            @inbounds out3[p, j, q] = β * A3[p, i, q] + α * A3[p, i + 1, q]
        end
    end
    return out3
end

# General fallback for external interpolators (DataInterpolations.jl etc.): one scalar
# interpolator per channel over a view, so nothing is copied and any `interp` works
# (vector-valued splines need slice types whose arithmetic round-trips).
function _tinterp_nd(interp, A, old_times, new_times, d; kws...)
    A3 = _as3d(A, d)
    fs = [interp(view(A3, p, :, q), old_times; kws...) for q in axes(A3, 3) for p in axes(A3, 1)]
    return _stack_channels(A, d, [f.(new_times) for f in fs])
end

# Assemble per-channel results (ordered as `vec` of the non-time dims) into `A`'s layout.
function _stack_channels(A, d, cols)
    out = similar(A, eltype(first(cols)), ntuple(i -> i == d ? length(first(cols)) : size(A, i), ndims(A)))
    out3 = _as3d(out, d)
    for (k, I) in enumerate(CartesianIndices((axes(out3, 1), axes(out3, 3))))
        out3[I[1], :, I[2]] = cols[k]
    end
    return out
end

function tinterp(A, t; dim = nothing, kws...)
    d = dimnum(A, dim)
    out = tinterp(unwrap(A), axiskeys(A, d), t; dim = d, kws...)
    return t isa AbstractArray ? rebuild(A, out, d, t) : out
end

"""
    tresample(A, old_times, freq; kw...)

Resample time series `A` onto a regular time grid with the specified frequency `freq`.

See also: [`tinterp`](@ref), [`time_grid`](@ref)
"""
tresample(A, old_times, freq; kw...) = tinterp(A, old_times, time_grid(old_times, freq); kw...)

function tresample(A, dt; dim = nothing, kws...)
    d = dimnum(A, dim)
    return tinterp(A, time_grid(axiskeys(A, d), dt); dim = d, kws...)
end


"""
    tsync(A, Bs...)

Synchronize multiple time series to have the same time points.

This function aligns time series `Bs...` to match time points of `A` by:

 1. Finding common time range between all time series
 2. Extracting subset of `A` within common range
 3. Interpolating each series in `Bs...` to match the time points of the subset of `A`

# Examples

```julia
A_sync, B_sync, C_sync = tsync(A, B, C)
```

See also: [`tinterp`](@ref), [`common_timerange`](@ref)
"""
function tsync(A, Bs...)
    tr = common_timerange(A, Bs...)
    @assert !isnothing(tr) "No common time range found"
    A_tsync = tclip(A, tr...)
    tstamps = times(A_tsync)
    Bs_syncs = map(Bs) do B
        tinterp(B, tstamps)
    end
    return A_tsync, Bs_syncs...
end


"""
    tinterp_nans(A; dim = nothing, kwargs...)

Interpolate only the NaN values in `A` along dimension `dim`.
"""
function tinterp_nans(A; dim = nothing, kwargs...)
    dims = dimnum(A, dim)
    t = axiskeys(A, dims)
    out = mapslices(parent(A); dims) do slice
        interpolate_nans!(slice, t; kwargs...)
    end
    return rebuild(A, out, dims, t)
end

# Interpolate only the NaN values in `u` along `t`.
function interpolate_nans!(u, t; interp = LinearInterpolation)
    # For 1D arrays, directly interpolate the NaN values
    nan_indices = findall(isnan, u)
    if !isempty(nan_indices) && length(nan_indices) < length(u)
        valid_indices = findall(!isnan, u)
        interp_obj = @views interp(u[valid_indices], t[valid_indices])
        for idx in nan_indices
            u[idx] = interp_obj(t[idx])
        end
    end
    return u
end


# Keeps `dt` disjoint from time vectors so the optional positional `dt` dispatches unambiguously.
const _Step = Union{Dates.Period, Real}

"""
    tfill_gaps(t, [dt]; max_gap = nothing) -> t_new
    tfill_gaps(A, [dt]; fill = NaN, max_gap = nothing, dim = nothing) -> typeof(A)
    tfill_gaps(A, t, [dt]; fill = NaN, max_gap = nothing, dim = ndims(A)) -> (A_new, t_new)

Insert missing timestamps (new samples set to `fill`); originals are kept exactly.
A gap `Δ` gets `round(Δ / dt) - 1` points at `t[i] + k * dt`; gaps `Δ > max_gap` are left open.
`dt` defaults to [`cadence`](@ref). Float eltypes are kept; others promote to fit `fill`.

```julia
tfill_gaps(da, Second(1); max_gap = Minute(5)) |> tinterp_nans   # fill short gaps, then interpolate
```
"""
function tfill_gaps(t::AbstractVector{<:AbstractTime}, dt::_Step = cadence(t; check = false); max_gap = nothing)
    fillers, n_new = _gap_fillers(t, dt, max_gap)
    return _fill_timestamps(t, fillers, n_new, dt)
end

function tfill_gaps(A, dt::Union{_Step, Nothing} = nothing; dim = nothing, kws...)
    d = dimnum(A, dim)
    out, t_new = _tfill_gaps(unwrap(A), axiskeys(A, d), d, dt; kws...)
    return rebuild(A, out, d, t_new)
end

tfill_gaps(A::AbstractArray, t::AbstractVector, dt::Union{_Step, Nothing} = nothing; dim = ndims(A), kws...) =
    _tfill_gaps(A, t, dim, dt; kws...)

# Counting first and writing into exactly-sized buffers beats a single push!-pass.
function _gap_fillers(t, dt, max_gap)
    dt_s = _seconds(dt)
    dt_s > 0 || throw(ArgumentError("dt must be positive, got $dt"))
    max_gap_s = isnothing(max_gap) ? Inf : _seconds(max_gap)
    n_old = length(t)
    fillers = Vector{Int}(undef, max(n_old - 1, 0))
    n_new = n_old
    @inbounds for i in 1:(n_old - 1)
        gap_s = _seconds(t[i + 1] - t[i])
        n = gap_s > max_gap_s ? 0 : max(0, round(Int, gap_s / dt_s) - 1)
        fillers[i] = n
        n_new += n
    end
    return fillers, n_new
end

function _fill_timestamps(t_old, fillers, n_new, dt)
    n_new == length(t_old) && return collect(t_old)
    t_new = Vector{eltype(t_old)}(undef, n_new)
    j = 1
    @inbounds for i in eachindex(fillers)
        ti = t_old[i]
        t_new[j] = ti
        for k in 1:fillers[i]
            t_new[j + k] = ti + k * dt
        end
        j += fillers[i] + 1
    end
    @inbounds t_new[end] = t_old[end]
    return t_new
end

_fill_eltype(::Type{T}, fill) where {T} = T <: AbstractFloat && fill isa Real ? T : promote_type(T, typeof(fill))

function _tfill_gaps(A::AbstractArray, t_old, d, dt; fill = NaN, max_gap = nothing)
    dt = isnothing(dt) ? cadence(t_old; check = false) : dt
    size(A, d) == length(t_old) || throw(DimensionMismatch("size(A, $d) = $(size(A, d)) but length(t) = $(length(t_old))"))
    fillers, n_new = _gap_fillers(t_old, dt, max_gap)
    t_new = _fill_timestamps(t_old, fillers, n_new, dt)
    sz = ntuple(i -> i == d ? n_new : size(A, i), ndims(A))
    out = similar(A, _fill_eltype(eltype(A), fill), sz)
    Base.fill!(out, fill)
    _copy_spread!(_as3d(out, d), _as3d(A, d), fillers)
    return out, t_new
end

# Copy slice `i` of `A3` to slice `i + sum(fillers[1:i-1])` of `out3`.
function _copy_spread!(out3, A3, fillers)
    offset = 0
    for i in axes(A3, 2)
        j = i + offset
        _foreach_pq(A3) do p, q
            @inbounds out3[p, j, q] = A3[p, i, q]
        end
        i <= length(fillers) && (offset += fillers[i])
    end
    return out3
end
