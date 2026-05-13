"""
    tfilter(A, times, Wn1, Wn2=nothing; dim=ndims(A), designmethod=nothing)
    tfilter(A, Wn1, Wn2=nothing; dim=nothing, designmethod=nothing)

Bandpass filter `A` between `Wn1` and `Wn2`. The upper cutoff defaults to
the Nyquist frequency.

References
- https://docs.juliadsp.org/stable/filters/
- https://www.mathworks.com/help/signal/ref/filtfilt.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
"""
function tfilter end

function tfilter(A, Wn1, Wn2 = nothing; dim = nothing, kwargs...)
    d = dimnum(A, dim)
    ts = axiskeys(A, d)
    data = tfilter(unwrap(A), ts, Wn1, Wn2; dim = d, kwargs...)
    return rebuild(A, data, d, ts)
end

"""
    pspectrum(A, times; dim=ndims(A), nfft=256, noverlap=div(nfft, 2), window=DSP.hamming)
    pspectrum(A; dim=nothing, freqdim=:frequency, kwargs...)

Compute a time-frequency power spectrum with `DSP.spectrogram`.

For containers that implement the common axis API, `times` are inferred from
the selected axis and the result is rebuilt with frequency and spectrogram time
axes followed by the remaining axes.
"""
function pspectrum end

function pspectrum(A; dim = nothing, freqdim = :frequency, kwargs...)
    d = dimnum(A, dim)
    ts = axiskeys(A, d)
    spec = pspectrum(unwrap(A), ts; dim = d, kwargs...)
    odims = Tuple(i for i in 1:ndims(A) if i != d)
    newdims = (freqdim, dims(A, d), ntuple(i -> dims(A, odims[i]), length(odims))...)
    keys = (spec.freq, spec.time, ntuple(i -> axiskeys(A, odims[i]), length(odims))...)
    return rebuild(A, spec.power, newdims, keys)
end
