"""
    dropna(A; dim=nothing)

Remove slices containing NaN values along dimension `dim`.
"""
@inline function dropna(A; dim = nothing)
    d = dimnum(A, dim)
    idxs = vec(all(!isnan, unwrap(A); dims = other_dims(A, d)))
    return @view A[_selectors(A, d, idxs)...]
end

"""
    tinterp_nans(A; dim = nothing, kwargs...)

Interpolate only the NaN values in `A` along dimension `dim`.
"""
tinterp_nans(A; dim = nothing, kwargs...) = _ontime(A, dim) do data, t, d
    x = rawview(t)
    out = mapslices(data; dims = d) do slice
        interpolate_nans!(slice, x; kwargs...)
    end
    return out, t
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


"""
    tfill_gaps(t, [dt]; max_gap = nothing) -> t_new
    tfill_gaps(A, [dt]; fill = NaN, max_gap = nothing, dim = nothing)

Insert missing timestamps (new samples set to `fill`); originals are kept exactly.
A gap `Δ` gets `round(Δ / dt) - 1` points at `t[i] + k * dt`; gaps `Δ > max_gap` are left open.
`dt` defaults to [`cadence`](@ref). Float eltypes are kept; others promote to fit `fill`.

```julia
tfill_gaps(da, Second(1); max_gap = Minute(5)) |> tinterp_nans   # fill short gaps, then interpolate
```
"""
function tfill_gaps(t::AbstractVector{<:AbstractTime}, dt = cadence(t; check = false); max_gap = nothing)
    fillers, n_new = _gap_fillers(t, dt, max_gap)
    return _fill_timestamps(t, fillers, n_new, dt)
end

tfill_gaps(A, dt = nothing; dim = nothing, kws...) = _ontime((data, t, d) -> _tfill_gaps(data, t, d, dt; kws...), A, dim)

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
