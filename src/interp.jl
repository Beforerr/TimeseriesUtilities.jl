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

struct LinearInterpolation end

function check_interp_times(t)
    issorted(t) || throw(ArgumentError("t must be sorted"))
    isempty(t) && throw(ArgumentError("at least one interpolation point is required"))
    return nothing
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

"""
    tinterp(A, t; interp=LinearInterpolation, dim=nothing, kws...)

Interpolate `A` at time `t`, dropping dimension `dim`, or at each of times `t`, keyed by `t`.

Any other `interp(u, t; kws...)` must return a callable; `kws` go to it.
"""
tinterp(A, t; kws...) = _tinterp(rebuild, A, t; kws...)

function _tinterp(rebuild::F, A, t; interp = LinearInterpolation, dim = nothing, kws...) where {F}
    # not `_ontime`: a closure would capture `interp` (a type, by default) as `DataType`, so the `_tinterp_nd` call would not infer
    d = resolve_dim(A, dim)
    data, t0 = unwrap(A), times(A, d)
    t isa AbstractArray && return rebuild(A, _interpolate(interp, data, t0, t, d; kws...), d, t)
    out = _interpolate(interp, data, t0, [t], d; kws...)
    return ndims(out) == 1 ? only(out) : _dropdim(A, dropdims(out; dims = d), d)
end

_interpolate(interp, A, t, tq, d; max_gap = nothing, kws...) =
    isnothing(max_gap) ? _tinterp_nd(interp, A, t, tq, d; kws...) : _tinterp_gapped(interp, A, t, tq, d, max_gap; kws...)

function _tinterp_gapped(interp, A, t, tq, d, max_gap; kws...)
    t, tq = _commontimes(t, tq)
    extrapolating = get(kws, :extrapolation, false) !== false
    gaps = Int[]
    hint = 0
    for (j, x) in enumerate(tq)
        i = hint = _interp_segment(t, x, true, hint)
        _ingap(t, x, i, extrapolating, max_gap) && push!(gaps, j)
    end
    if length(gaps) == length(tq) # no interpolator: some reject empty or too few samples
        T = float(eltype(A))
        return fill!(similar(A, T, ntuple(i -> i == d ? length(tq) : size(A, i), ndims(A))), zero(T) / 0)
    end
    qs = collect(promote_type(eltype(t), eltype(tq)), tq)
    qs[gaps] .= clamp.(view(qs, gaps), first(t), last(t))
    out = float(_tinterp_nd(interp, A, t, qs, d; kws...))
    _as3d(out, d)[:, gaps, :] .= zero(eltype(out)) / 0
    return out
end

function _ingap(t, x, i, extrapolating, max_gap)
    isempty(t) && return true
    x < first(t) && return !extrapolating || first(t) - x > max_gap
    x > last(t) && return !extrapolating || x - last(t) > max_gap
    return length(t) > 1 && x != t[i] && x != t[i + 1] && t[i + 1] - t[i] > max_gap
end

# Fast one pass over all channels per query point
@inline _tinterp_nd(::Type{LinearInterpolation}, A, old_times, new_times, d; extrapolation = false, check = true) =
    _tinterp_linear_nd(A, old_times, new_times, d, extrapolation, check)

function _tinterp_linear_nd(A, old_times, new_times, d, extrapolation, check)
    old_times, new_times = _commontimes(old_times, new_times)
    out_sz = ntuple(i -> i == d ? length(new_times) : size(A, i), ndims(A))
    out = similar(A, float(eltype(A)), out_sz)
    old_times == new_times && return copyto!(out, A)
    check && check_interp_times(old_times)
    _tinterp_linear!(_samples(out, d), _samples(A, d), old_times, new_times, extrapolation)
    return out
end

function _tinterp_linear!(out, A, old_times, new_times, extrapolation)
    hint = 0
    @inbounds for (j, x) in enumerate(new_times)
        i = hint = _interp_segment(old_times, x, extrapolation, hint)
        k = length(old_times) == 1 || x == old_times[i] ? i : x == old_times[i + 1] ? i + 1 : 0
        if k > 0
            _foreach_channel(A) do p, q
                @inbounds out[p..., j, q...] = A[p..., k, q...]
            end
            continue
        end
        α = (x - old_times[i]) / (old_times[i + 1] - old_times[i])
        β = 1 - α
        _foreach_channel(A) do p, q
            @inbounds out[p..., j, q...] = β * A[p..., i, q...] + α * A[p..., i + 1, q...]
        end
    end
    return out
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

Interpolate `A` onto `tmin:dt:tmax` of its times.
"""
function tresample(A, dt; dim = nothing, kws...)
    return _tinterp(_rebuild_ordered, A, _time_grid(times(A, dim), dt); dim, kws...)
end


"""
    tsync(A, Bs...; dim=nothing, kws...)

Clip `A` to the time range common to all series and interpolate each of `Bs` onto its times.
Throws `ArgumentError` if the ranges do not overlap.
"""
function tsync(A, Bs...; dim = nothing, kws...)
    tr = @something trange(A, Bs...; dim) throw(ArgumentError("No common time range found"))
    A_tsync = tclip(A, tr; dim)
    tstamps = times(A_tsync, dim)
    Bs_syncs = map(Bs) do B
        _tinterp(_rebuild_ordered, B, tstamps; dim, kws...)
    end
    return A_tsync, Bs_syncs...
end

