module TimeseriesUtilitiesDSPExt

using DSP
using DSP: Butterworth, Bandpass, digitalfilter, filtfilt
import TimeseriesUtilities as TU
using TimeseriesUtilities: TimeOffsets, samplingrate

function TU._tfilter(A::AbstractArray, ts::AbstractVector, d::Integer, Wn1::Real, Wn2; designmethod = Butterworth(2), max_gap = 0)
    dt = TU.cadence(ts)
    fs = inv(TU._seconds(dt))
    Wn2 = @something(Wn2, 0.999 * fs / 2)
    f = digitalfilter(Bandpass(Wn1, Wn2), designmethod; fs)
    limit = isnothing(max_gap) ? Inf : TU._seconds(max_gap)
    limit >= 0 || throw(ArgumentError("max_gap must be nonnegative"))
    # Track missing samples without allocating across long gaps.
    grid = _sample_grid(ts, TU._seconds(dt))
    return mapslices(A; dims = d) do x
        last(grid) == length(x) && !any(isnan, x) && return filtfilt(f, x)
        _filter_gaps(f, x, ts, grid, limit)
    end
end

_sample_grid(ts::AbstractRange, dt_s) = eachindex(ts)
function _sample_grid(ts, dt_s)
    grid = Vector{Int}(undef, length(ts))
    grid[1] = 1
    for i in 2:length(ts)
        gap_s = TU._seconds(ts[i] - ts[i - 1])
        gap_s > 0 || throw(ArgumentError("filter timestamps must increase"))
        grid[i] = grid[i - 1] + 1 + TU._missing_steps(gap_s, dt_s)
    end
    return grid
end

function _filter_gaps(f, x, ts, grid, limit)
    valid = findall(!isnan, x)
    isempty(valid) && return filtfilt(f, x)
    out = nothing
    start = 1
    for k in eachindex(valid)
        if k < length(valid)
            i, j = valid[k], valid[k + 1]
            if grid[j] == grid[i] + 1 || TU._seconds(ts[j] - ts[i]) <= limit
                continue
            end
        end
        indices = @view valid[start:k]
        filtered = _filter_run(f, x, grid, indices)
        isnothing(out) && (out = fill!(similar(x, eltype(filtered)), NaN))
        out[indices] = filtered
        start = k + 1
    end
    return out
end

function _filter_run(f, x, grid, indices)
    first_slot = grid[first(indices)]
    n = grid[last(indices)] - first_slot + 1
    n == length(indices) && return filtfilt(f, view(x, first(indices):last(indices)))
    slots = view(grid, indices) .- (first_slot - 1)
    u = fill!(similar(x, float(eltype(x)), n), NaN)
    u[slots] = view(x, indices)
    TU.interpolate_nans!(u, eachindex(u))
    return filtfilt(f, u)[slots]
end

function _spectrogram(x, ts, fs; nfft = 256, noverlap = div(nfft, 2), window = DSP.hamming)
    spec = DSP.spectrogram(x, nfft, noverlap; fs, window)
    return (power = spec.power, time = TimeOffsets(spec.time, first(ts)), freq = spec.freq)
end

function TU._pspectrum(A::AbstractArray, ts::AbstractVector, d::Integer; kw...)
    fs = samplingrate(ts)
    ndims(A) == 1 && return _spectrogram(A, ts, fs; kw...)
    specs = map(slice -> _spectrogram(vec(slice), ts, fs; kw...), eachslice(A; dims = TU.other_dims(A, d)))
    first_spec = first(specs)
    return (power = stack(spec.power for spec in specs), time = first_spec.time, freq = first_spec.freq)
end

end
