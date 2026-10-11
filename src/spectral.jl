"""
    tfilter(A, Wn1, Wn2=nothing; dim=nothing, designmethod=Butterworth(2))

Zero-phase bandpass filter `A` between `Wn1` and `Wn2` hertz; the upper cutoff defaults to just below Nyquist.
Requires DSP.jl.

References
- https://docs.juliadsp.org/stable/filters/
- https://www.mathworks.com/help/signal/ref/filtfilt.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
"""
tfilter(A, Wn1, Wn2 = nothing; dim = nothing, kw...) =
    _ontime((data, t, d) -> (_tfilter(data, t, d, Wn1, Wn2; kw...), t), A, dim)

"""
    pspectrum(A; dim=nothing, freqdim=:frequency, nfft=256, noverlap=div(nfft, 2), window=DSP.hamming)

Short-time Fourier power spectrum of `A` with `DSP.spectrogram`, with axes
(frequency, spectrogram time, remaining axes of `A`). Requires DSP.jl.
"""
function pspectrum(A; dim = nothing, freqdim = :frequency, kwargs...)
    d = dimnum(A, dim)
    spec = _pspectrum(unwrap(A), axiskeys(A, d), d; kwargs...)
    odims = other_dims(A, d)
    newdims = (freqdim, Val(d), map(Val, odims)...)
    keys = (spec.freq, spec.time, map(i -> axiskeys(A, i), odims)...)
    return rebuild(A, spec.power, newdims, keys)
end

# Kernels on plain data `A` with times `t` along dimension `d`, defined by the DSP extension.
_tfilter(A, t, d, Wn1, Wn2; kw...) = throw(ArgumentError("`tfilter` requires `using DSP`"))
_pspectrum(A, t, d; kw...) = throw(ArgumentError("`pspectrum` requires `using DSP`"))
