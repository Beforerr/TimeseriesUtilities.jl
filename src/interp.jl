rawview(x) = x
rawview(x::AbstractTime) = Dates.value(x)
rawview(x::AbstractArray{<:AbstractTime}) = Dates.value.(x)
# Zero-copy, and indexes as fast as the broadcast copy; assumes `Dates.value` returns the innermost field.
function rawview(x::DenseArray{T}) where {T <: AbstractTime}
    R = Base.promote_op(Dates.value, T)
    return _wraps(T, R) ? reinterpret(R, x) : Dates.value.(x)
end
_wraps(T, R) = T === R || (isbitstype(T) && fieldcount(T) == 1 && _wraps(fieldtype(T, 1), R))

# Mixed time types (`Date` knots, `DateTime` queries) neither subtract nor share raw-value
# epochs and units, so both sides convert to the promoted time type first.
function _commontimes(t, tq)
    T = promote_type(eltype(t), eltype(tq))
    T <: AbstractTime && !isconcretetype(T) && eltype(t) != eltype(tq) &&
        throw(ArgumentError("times of types $(eltype(t)) and $(eltype(tq)) have no common type"))
    return _astime(T, t), _astime(T, tq)
end
_astime(::Type, x) = x
_astime(::Type{T}, x) where {T <: AbstractTime} = eltype(x) <: T ? x : convert.(T, x)

# Numeric interpolators (DataInterpolations.jl, FastInterpolations.jl) cannot take `Dates` times
_rawtimes(t, tq) = map(rawview, _commontimes(t, tq))

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
    tinterp(A, t; interp=LinearInterpolation, dim=nothing, kws...)

Interpolate `A` at time `t`, or at each of times `t` (keyed by `t`).

The `interp` constructor must accept `(u, t; kws...)` and return a callable object.

# Examples

```julia
# Interpolate at a single time point
tinterp(time_series, DateTime("2023-01-01T12:00:00"))

# Interpolate at multiple time points using cubic spline interpolation
new_times = DateTime("2023-01-01"):Hour(1):DateTime("2023-01-02")
tinterp(time_series, new_times; interp = CubicSpline)
```
"""
function tinterp(A, t; interp = LinearInterpolation, dim = nothing, kws...)
    # not `_ontime`: a closure would capture `interp` (a type, by default) as `UnionAll`, so the `_tinterp_nd` call would not infer
    d = resolve_dim(A, dim)
    data, t0 = unwrap(A), axiskeys(A, d)
    t isa AbstractArray && return rebuild(A, _tinterp_nd(interp, data, t0, t, d; kws...), d, t)
    out = _tinterp_nd(interp, data, t0, [t], d; kws...)
    return ndims(out) == 1 ? only(out) : _dropdim(A, dropdims(out; dims = d), d)
end

# Fast path for built-in LinearInterpolation: one pass over all channels per query point.
@inline _tinterp_nd(::Type{<:LinearInterpolation}, A, old_times, new_times, d; extrapolation = false, kws...) =
    _tinterp_linear_nd(A, old_times, new_times, d, extrapolation)

function _tinterp_linear_nd(A, old_times, new_times, d, extrapolation)
    old_times, new_times = _commontimes(old_times, new_times)
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
    x, xq = _rawtimes(old_times, new_times)
    A3 = _as3d(A, d)
    ctor, kw = _callable(interp), _interp_kws(interp; kws...)
    fs = [ctor(view(A3, p, :, q), x; kw...) for q in axes(A3, 3) for p in axes(A3, 1)]
    return _stack_channels(A, d, [f.(xq) for f in fs])
end

# Keywords for constructor `interp`; extensions map a `Bool` `extrapolation` to the constructor's own option.
_interp_kws(interp; kws...) = kws

# A closure capturing type values stores them as `DataType`/`UnionAll` fields
struct _Ctor{T} end
(::_Ctor{T})(args...; kws...) where {T} = T(args...; kws...)
_callable(T::Type) = _Ctor{T}()
_callable(f) = f

"""
    tresample(A, dt; dim=nothing, kws...)

Interpolate `A` onto the regular time grid of step `dt` spanning its times.

See also: [`tinterp`](@ref)
"""
function tresample(A, dt; dim = nothing, kws...)
    return tinterp(A, _time_grid(times(A, dim), dt); dim, kws...)
end


"""
    tsync(A, Bs...; dim=nothing, kws...)

Synchronize multiple time series to have the same time points.

This function aligns time series `Bs...` to match time points of `A` by:

 1. Finding common time range between all time series
 2. Extracting subset of `A` within common range
 3. Interpolating each series in `Bs...` to match the time points of the subset of `A`, with [`tinterp`](@ref)`(B, t; kws...)`

# Examples

```julia
A_sync, B_sync, C_sync = tsync(A, B, C)
```
"""
function tsync(A, Bs...; dim = nothing, kws...)
    tr = @something trange(A, Bs...; dim) throw(ArgumentError("No common time range found"))
    A_tsync = tclip(A, tr; dim)
    tstamps = times(A_tsync, dim)
    Bs_syncs = map(Bs) do B
        tinterp(B, tstamps; dim, kws...)
    end
    return A_tsync, Bs_syncs...
end

