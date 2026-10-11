module TimeseriesUtilitiesDSPExt

using DSP
using DSP: Butterworth, Bandpass, digitalfilter, filtfilt
import TimeseriesUtilities as TU
using TimeseriesUtilities: TimeOffsets, samplingrate

function TU._tfilter(A::AbstractArray, ts::AbstractVector, d::Integer, Wn1::Real, Wn2; designmethod = Butterworth(2))
    fs = samplingrate(ts)
    Wn2 = @something(Wn2, 0.999 * fs / 2)
    f = digitalfilter(Bandpass(Wn1, Wn2), designmethod; fs)
    return mapslices(slice -> filtfilt(f, slice), A; dims = d)
end

function _spectrogram(x, ts; nfft = 256, noverlap = div(nfft, 2), window = DSP.hamming)
    spec = DSP.spectrogram(x, nfft, noverlap; fs = samplingrate(ts), window)
    return (power = spec.power, time = TimeOffsets(spec.time, first(ts)), freq = spec.freq)
end

function TU._pspectrum(A::AbstractArray, ts::AbstractVector, d::Integer; kw...)
    ndims(A) == 1 && return _spectrogram(A, ts; kw...)
    specs = map(slice -> _spectrogram(vec(slice), ts; kw...), eachslice(A; dims = TU.other_dims(A, d)))
    first_spec = first(specs)
    return (power = stack(spec.power for spec in specs), time = first_spec.time, freq = first_spec.freq)
end

end
