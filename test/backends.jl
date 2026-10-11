@testitem "backends agree" begin
    using AxisKeys: KeyedArray
    using DimensionalData: DimArray, Ti, Y
    using DataInterpolations: QuadraticInterpolation
    using TimeseriesUtilities: unwrap

    A = [1.0 10.0; NaN 20.0; 3.0 30.0; 4.0 400.0; 5.0 50.0]
    t = [0.0, 1.0, 2.0, 4.0, 5.0]
    backends = (
        KeyedArray = KeyedArray(A; time = t, c = 1:2),
        DimArray = DimArray(A, (Ti(t), Y(1:2))),
    )
    ops = (
        x -> tclip(x, 1.0, 4.0),
        x -> tview(x, 1.0, 4.0),
        x -> tselect(x, 1.2),
        x -> tmask(x, 1.0, 2.0),
        x -> tshift(x, 1.0),
        dropna,
        x -> smooth(x, 2.0),
        tderiv,
        tsubtract,
        x -> tinterp(x, [0.5, 4.5]),
        x -> tinterp(x, [0.5, 4.5]; interp = QuadraticInterpolation),
        x -> tresample(x, 0.5),
        x -> tfill_gaps(x, 1.0),
        tinterp_nans,
        tmean,
        x -> tmean(x, 2.0; period = 3.0),
        x -> mapsamples(sum, x),
        x -> outlier_mask(x, 3.0; threshold = 1),
        x -> despike(x; threshold = 1),
    )
    # 1.10 loses inference in several kernels and in DimensionalData
    run(op, x) = VERSION >= v"1.11" ? @inferred(op(x)) : op(x)
    @testset "op $i on $name" for (i, op) in enumerate(ops), (name, x) in pairs(backends)
        ref = run(op, Timed(A, t; dim = 1))
        r = run(op, x)
        @test isequal(unwrap(r), unwrap(ref))
        ref isa Timed && @test times(r) == times(ref)
    end
end
