module TimeseriesUtilitiesUnitfulExt

using Dates
import TimeseriesUtilities as TU
using Unitful

function TU._time_grid(ts, dt::Unitful.Quantity)
    tmin, tmax = extrema(ts)
    return if dimension(dt) == Unitful.𝐓
        tmin:_2dates(dt):tmax
    elseif dimension(dt) == Unitful.𝐓^-1
        _dt = round(Nanosecond, 1 / dt)
        tmin:_dt:tmax
    else
        tmin:dt:tmax
    end
end

for (period, unit) in (
        (Dates.Week, Unitful.wk), (Dates.Day, Unitful.d), (Dates.Hour, Unitful.hr),
        (Dates.Minute, Unitful.minute), (Dates.Second, Unitful.s), (Dates.Millisecond, Unitful.ms),
        (Dates.Microsecond, Unitful.μs), (Dates.Nanosecond, Unitful.ns),
    )
    @eval _2dates(::typeof($unit)) = $period
end

_2dates(x::Unitful.Quantity) = _2dates(Unitful.unit(x))(x)

end
