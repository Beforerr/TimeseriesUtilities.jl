"""
    tfilter(da, Wn1, Wn2=nothing; designmethod=nothing)

Bandpass filter `da` between `Wn1` and `Wn2`. By default, the upper cutoff is the Nyquist frequency.

References
- https://docs.juliadsp.org/stable/filters/
- https://www.mathworks.com/help/signal/ref/filtfilt.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html

Issues
- DSP.jl and Unitful.jl: https://github.com/JuliaDSP/DSP.jl/issues/431
"""
function tfilter(da::AbstractDimArray, Wn1, Wn2 = nothing; designmethod = nothing)
    designmethod = @something(designmethod, Butterworth(2))
    fs = samplingrate(da)
    Wn2 = @something(Wn2, 0.999 * fs / 2)
    Wn1, Wn2, fs = (Wn1, Wn2, fs) ./ 1u"Hz" .|> NoUnits
    f = digitalfilter(Bandpass(Wn1, Wn2), designmethod; fs)
    res = filtfilt(f, ustrip(parent(da)))
    return rebuild(da; data = res * (da |> eltype |> unit))
end
