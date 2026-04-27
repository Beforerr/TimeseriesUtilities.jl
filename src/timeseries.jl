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
    smooth(data, times, window; dim=ndims(data))
    smooth(data, window; dim=ndims(data))

Smooths a time series by computing a moving average over a sliding window.
Edge windows are truncated, so the output has the same size as the input.

Each window covers the half-open interval `[coord - before, coord + after)`.
A scalar `window` is interpreted as a coordinate span along the smoothed axis.

# Arguments
- `dim=ndims(data)`: Dimension along which to perform smoothing
- `op=nanmean`: Function used to aggregate each window
"""
function smooth(data, coords, window; dim = ndims(data), op = nanmean)
    length(coords) == size(data, dim) || throw(DimensionMismatch("length(coords) must match size(data, dim)"))
    issorted(coords) || throw(ArgumentError("coords must be sorted"))
    before, after = _window_offsets(window)
    return mapslices(data; dims = dim) do slice
        op.(SlidingWindow(slice, coords, before, after))
    end
end

function smooth(data, window; dim = ndims(data), kw...)
    return smooth(data, axes(data, dim), window; dim, kw...)
end


function _window_offsets(window::Tuple)
    @assert length(window) == 2
    return window[1], window[2]
end

_window_offsets(window) = _half(window), _half(window)
_half(window) = window / 2
_half(window::Period) = Millisecond(window) / 2

"""
    tfilter(da, Wn1, Wn2=samplingrate(da) / 2; designmethod=nothing)

By default, the max frequency corresponding to the Nyquist frequency is used.

References
- https://docs.juliadsp.org/stable/filters/
- https://www.mathworks.com/help/signal/ref/filtfilt.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html

Issues
- DSP.jl and Unitful.jl: https://github.com/JuliaDSP/DSP.jl/issues/431
"""
function tfilter(da::AbstractDimArray, Wn1, Wn2 = 0.999 * samplingrate(da) / 2; designmethod = nothing)
    designmethod = something(designmethod, Butterworth(2))
    fs = samplingrate(da)
    Wn1, Wn2, fs = (Wn1, Wn2, fs) ./ 1u"Hz" .|> NoUnits
    f = digitalfilter(Bandpass(Wn1, Wn2; fs), designmethod)
    res = filtfilt(f, ustrip(parent(da)))
    return rebuild(da; data = res * (da |> eltype |> unit))
end

function _dropna(A; dim = nothing)
    valid_idx = vec(all(!isnan, A; dims = other_dims(A, dim)))
    return selectdim(A, dim, valid_idx)
end

"""
    dropna(A; dim=nothing)
    dropna(A::AbstractDimArray; dim=nothing, query=nothing)

Remove slices containing NaN values along along the `dim` dimension.
"""
dropna(A; dim = nothing) = _dropna(A; dim)

function dropna(A::AbstractDimArray; query = nothing, dim = nothing)
    dim = @something dim dimnum(A, query)
    return _dropna(A; dim)
end

function dropna(ds::DimStack, query = nothing)
    query = something(query, TimeDim)
    Dim, T = dimtype_eltype(ds, query)
    dims = otherdims(ds, query)

    valid_idx = mapreduce(.*, values(ds)) do A
        vec(all(!isnan, A; dims))
    end

    return ds[Dim(valid_idx)]
end


function rectify(ts::DimensionalData.Dimension; tol = 4, atol = nothing)
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

"""Rectify the time step of a `DimArray` to be uniform."""
function rectify(da; tol = 2, kwargs...)
    times = dims(da, Ti)
    t0 = times[1]
    dtime = Quantity.(times.val .- t0)
    new_times = rectify(Ti(dtime); tol)
    return set(da, Ti => new_times .+ t0)
end

# """
#     tsplit(da::AbstractDimArray, dim=Ti)

# Splits up data along dimension `dim`.
# """
# function tsplit(da::AbstractDimArray, dim = Ti; new_names = labels(da))
#     odims = otherdims(da, dim)
#     rows = eachslice(da; dims = odims)
#     das = map(rows, new_names) do row, name
#         rename(modify_meta(row; long_name = name), name)
#     end
#     return DimStack(das...)
# end

for f in (:smooth, :tfilter)
    @eval $f(args...; kwargs...) = da -> $f(da, args...; kwargs...)
end
