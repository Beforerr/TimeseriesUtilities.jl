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
@inline function tinterp(A, old_times, new_times; interp = LinearInterpolation, dim = ndims(A), kws...)
    return if ndims(A) == 1
        interp(A, old_times; kws...).(new_times)
    else
        f = interp(LazyStaticSlices(A, dim), old_times; kws...)
        stack(f, new_times; dims = dim)
    end
end

function tinterp(A, t; dim = nothing, kws...)
    d = dimnum(A, dim)
    out = tinterp(unwrap(A), axiskeys(A, d), t; dim = d, kws...)
    return t isa AbstractArray ? rebuild_axis(A, out, d, t) : out
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
    return rebuild_axis(A, out, dims, t)
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
