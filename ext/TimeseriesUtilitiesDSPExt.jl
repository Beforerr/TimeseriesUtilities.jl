module TimeseriesUtilitiesDSPExt

using DSP
import TimeseriesUtilities
using TimeseriesUtilities: samplingrate
using DimensionalData: AbstractDimArray, rebuild
using Unitful

function TimeseriesUtilities.tfilter(da::AbstractDimArray, Wn1, Wn2 = nothing; designmethod = nothing)
    designmethod = @something(designmethod, Butterworth(2))
    fs = samplingrate(da)
    Wn2 = @something(Wn2, 0.999 * fs / 2)
    Wn1, Wn2, fs = (Wn1, Wn2, fs) ./ 1u"Hz" .|> NoUnits
    f = digitalfilter(Bandpass(Wn1, Wn2), designmethod; fs)
    res = filtfilt(f, ustrip.(parent(da)))
    return rebuild(da; data = res * (da |> eltype |> unit))
end

end
