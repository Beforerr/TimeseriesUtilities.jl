"""
    tfilter(A, Wn1, Wn2=nothing; dim=nothing, designmethod=Butterworth(2), max_gap=0)

Apply a zero-phase bandpass filter between `Wn1` and `Wn2` hertz using DSP.jl.
The upper cutoff defaults to just below Nyquist.

NaNs and omitted timestamps split the signal into runs filtered independently.
To bridge a gap by linear interpolation, set `max_gap` to the maximum time
between bracketing observations (seconds or a time period), or `nothing` for
no limit. The output retains the original timestamps and NaNs.
"""
tfilter(A, Wn1, Wn2 = nothing; dim = nothing, kw...) =
    _ontime((data, t, d) -> (_tfilter(data, t, d, Wn1, Wn2; kw...), t), A, dim)

"""
    pspectrum(A; dim=nothing, freqdim=:frequency, nfft=256, noverlap=div(nfft, 2), window=DSP.hamming)

Short-time Fourier power spectrum of `A` with `DSP.spectrogram`, with axes
(frequency, spectrogram time, remaining axes of `A`). Requires DSP.jl.
"""
function pspectrum(A; dim = nothing, freqdim = :frequency, kwargs...)
    d = resolve_dim(A, dim)
    spec = _pspectrum(unwrap(A), times(A, d), d; kwargs...)
    odims = other_dims(A, d)
    return rebuild(A, spec.power, (freqdim => spec.freq, Val(d) => spec.time, map(Val, odims)...))
end

# Kernels on plain data `A` with times `t` along dimension `d`, defined by the DSP extension.
_tfilter(A, t, d, Wn1, Wn2; kw...) = throw(ArgumentError("`tfilter` requires `using DSP`"))
_pspectrum(A, t, d; kw...) = throw(ArgumentError("`pspectrum` requires `using DSP`"))
