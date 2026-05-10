module TimeseriesUtilitiesDSPExt

using DSP
using DSP: Butterworth, Bandpass, digitalfilter, filtfilt
import TimeseriesUtilities
using TimeseriesUtilities: TimeOffsets

function _sampling_frequency(ts)
    dt = TimeseriesUtilities.resolution(ts)
    return inv(TimeseriesUtilities._seconds(dt))
end

function TimeseriesUtilities.tfilter(
        A::AbstractArray,
        ts::AbstractVector,
        Wn1::Real,
        Wn2 = nothing;
        designmethod = Butterworth(2),
        dim = ndims(A),
    )
    length(ts) == size(A, dim) || throw(DimensionMismatch("length(times) must match size(data, dim)"))
    fs = _sampling_frequency(ts)
    Wn2 = @something(Wn2, 0.999 * fs / 2)
    f = digitalfilter(Bandpass(Wn1, Wn2), designmethod; fs)
    return mapslices(slice -> filtfilt(f, slice), A; dims = dim)
end

"""
    pspectrum(data, times; nfft=256, noverlap=div(nfft, 2), window=DSP.hamming)

Compute a short-time Fourier power spectrum from vector data and coordinates.
"""
function TimeseriesUtilities.pspectrum(
        x::AbstractVector,
        ts::AbstractVector;
        dim = nothing,
        nfft = 256,
        noverlap = div(nfft, 2),
        window = DSP.hamming,
    )
    dim in (nothing, 1) || throw(ArgumentError("dimension $dim out of range for vector input"))
    length(x) == length(ts) || throw(DimensionMismatch("data and times must have the same length"))
    spec = DSP.spectrogram(x, nfft, noverlap; fs = _sampling_frequency(ts), window)
    return (power = spec.power, time = TimeOffsets(spec.time, first(ts)), freq = spec.freq)
end

function TimeseriesUtilities.pspectrum(
        A::AbstractArray,
        ts::AbstractVector;
        dim = ndims(A),
        nfft = 256,
        noverlap = div(nfft, 2),
        window = DSP.hamming,
    )
    length(ts) == size(A, dim) || throw(DimensionMismatch("length(times) must match size(data, dim)"))
    ndims(A) == 1 && return TimeseriesUtilities.pspectrum(vec(A), ts; nfft, noverlap, window)

    odims = Tuple(i for i in 1:ndims(A) if i != dim)
    specs = map(eachslice(A; dims = odims)) do slice
        TimeseriesUtilities.pspectrum(vec(slice), ts; nfft, noverlap, window)
    end

    first_spec = first(specs)
    power = stack(spec.power for spec in specs)
    return (power = power, time = first_spec.time, freq = first_spec.freq)
end

end
