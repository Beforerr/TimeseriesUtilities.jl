module TimeseriesUtilitiesFastInterpolationsExt

using FastInterpolations
import TimeseriesUtilities as TU

# Families with a `Series` method: one interval search per query point, shared by all channels.
const SeriesInterp = Union{typeof(constant_interp), typeof(linear_interp), typeof(quadratic_interp), typeof(cubic_interp)}
const LocalInterp = Union{typeof(pchip_interp), typeof(akima_interp), typeof(cardinal_interp)}

_extrap(e::Bool) = e ? ExtendExtrap() : NoExtrap()
_extrap(e) = e

function TU._tinterp_nd(f::Union{SeriesInterp, LocalInterp}, A, old_times, new_times, d; extrapolation = false, kws...)
    x, xq = TU.rawview(old_times), TU.rawview(new_times)
    A3 = TU._as3d(A, d)
    kw = (; extrap = _extrap(extrapolation), kws...)
    cols = if f isa SeriesInterp
        Y = size(A3, 1) == 1 ? reshape(A3, size(A3, 2), :) : reshape(permutedims(A3, (2, 1, 3)), size(A3, 2), :)
        f(x, Series(Y), xq; kw...)
    else
        [f(x, view(A3, p, :, q), xq; kw...) for q in axes(A3, 3) for p in axes(A3, 1)]
    end
    return TU._stack_channels(A, d, cols)
end

end
