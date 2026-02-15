function resolution(times; tol = 2, check = true, f = _median!)
    dts = diff(times)
    return if check
        dt, relerr = stat_relerr(dts, f)
        relerr > exp10(-tol - 1) && @warn "Time resolution is is not approximately constant (relerr ≈ $relerr)"
        dt
    else
        f(dts)
    end
end

samplingrate(da) = 1u"s" / resolution(da) * u"Hz" |> u"Hz"


"""
    smooth(data, window::Integer; dim=1, kwargs...)
    smooth(data, times, window; dim=1, kwargs...)

Smooths a time series by computing a moving average over a sliding window.

The size of the sliding `window` can be either:
  - `Quantity`: A time duration that will be converted to number of samples based on data resolution
  - `Integer`: Number of samples directly

# Arguments
- `dim=1`: Dimension along which to perform smoothing
- `kwargs...`: Additional arguments passed to `RollingWindowArrays.rolling`
"""
function smooth(data::AbstractArray, window::Integer; dim = 1, kwargs...)
    return mapslices(data; dims = dim) do slice
        nanmean.(RollingWindowArrays.rolling(slice, window; kwargs...))
    end
end

smooth(data, times, window; kwargs...) = smooth(data, Integer(div(window, resolution(times))); kwargs...)

"""
    tfilter(data, fs, Wn1, Wn2=0.999*fs/2; designmethod=nothing)

Apply a bandpass filter to `data` with sampling rate `fs`.
By default, the max frequency corresponding to the Nyquist frequency is used.

References
- https://docs.juliadsp.org/stable/filters/
- https://www.mathworks.com/help/signal/ref/filtfilt.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html

Issues
- DSP.jl and Unitful.jl: https://github.com/JuliaDSP/DSP.jl/issues/431
"""
function tfilter(data::AbstractArray, fs, Wn1, Wn2 = 0.999 * fs / 2; designmethod = nothing)
    designmethod = something(designmethod, Butterworth(2))
    Wn1, Wn2, fs = (Wn1, Wn2, fs) ./ 1u"Hz" .|> NoUnits
    f = digitalfilter(Bandpass(Wn1, Wn2; fs), designmethod)
    res = filtfilt(f, ustrip(data))
    return res * (eltype(data) |> unit)
end

"""
    dropna(A; dim)

Remove slices containing NaN values along the `dim` dimension.
"""
function dropna(A; dim = nothing)
    valid_idx = vec(all(!isnan, A; dims = other_dims(A, dim)))
    return selectdim(A, dim, valid_idx)
end

"""
    rectify(ts::AbstractVector; tol=4, atol=nothing)

Rectify a time vector to have uniform step size.
"""
function rectify(ts::AbstractVector; tol = 4, atol = nothing)
    u = unit(eltype(ts))
    ts = collect(ts)
    stp = ts |> diff |> mean
    err = ts |> diff |> std
    tol = Int(tol - round(log10(stp |> ustripall)))

    if isnothing(atol) && ustripall(err) > exp10(-tol - 1)
        @warn "Step $stp is not approximately constant (err=$err, tol=$(exp10(-tol - 1))), skipping rectification"
    else
        if !isnothing(atol)
            tol = atol
        end
        stp = u == NoUnits ? round(stp; digits = tol) : round(u, stp; digits = tol)
        t0, t1 = u == NoUnits ? round.(extrema(ts); digits = tol) :
            round.(u, extrema(ts); digits = tol)
        ts = range(start = t0, step = stp, length = length(ts))
    end
    return ts
end

for f in (:smooth, :tfilter)
    @eval $f(args...; kwargs...) = da -> $f(da, args...; kwargs...)
end
