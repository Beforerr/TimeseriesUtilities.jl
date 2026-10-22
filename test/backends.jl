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
        x -> tsmooth(x, 2.0),
        tderiv,
        x -> tsubtract(first, x),
        x -> tsubtract(maximum, x, 2.0),
        x -> tinterp(x, [0.5, 4.5]),
        x -> tinterp(x, [0.5, 4.5]; interp = QuadraticInterpolation),
        x -> tresample(x, 0.5),
        x -> tfill_gaps(x, 1.0),
        tinterp_na,
        tmean,
        x -> tmean(x, 2.0; period = 3.0),
        x -> tmaximum(x, 2.0),
        targmax,
        x -> mapsamples(sum, x),
        x -> mapsamples(v -> reshape(v, 1, :), x),
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

@testitem "kept axes keep their lookup" begin
    using DimensionalData
    using DimensionalData.Lookups: Center, ForwardOrdered, Intervals, Irregular, Sampled, span

    e = Dim{:energy}(Sampled([1.0, 2.0, 5.0]; span = Irregular(0.5, 6.0), sampling = Intervals(Center()), order = ForwardOrdered()))
    x = DimArray(rand(4, 3), (Ti(1.0:4.0), e))
    @test span(dims(tmean(x), :energy)) == Irregular(0.5, 6.0)
end

@testitem "dimensions without coordinates" begin
    using AxisKeys: KeyedArray
    using DimensionalData: DimArray, Ti, Y, lookup, NoLookup
    using TimeseriesUtilities: axiskeys, rebuild, rekey, unwrap

    data = reshape(1:12, 4, 3)
    t = 10:10:40
    backends = (Timed(data, t; dim = 1), KeyedArray(data; time = t, c = 1:3), DimArray(data, (Ti(t), Y(1:3))))
    for x in backends
        # KeyedArray cannot represent an axis without coordinates, so it gets positions
        y = mapsamples(v -> v[1:2], x)
        @test isnothing(axiskeys(y, 2)) == !(x isa KeyedArray)
        @test isnothing(axiskeys(rebuild(y, unwrap(y)), 2)) == !(x isa KeyedArray)
        z = rekey(y, 2, [100, 200])
        @test axiskeys(z, 2) == [100, 200]
        @test_throws DimensionMismatch mapsamples(+, z, rekey(z, 2, [300, 400]); dim = 2)
        @test unwrap(mapsamples(+, x, data; dim = 1)) == 2data
        @test_throws DimensionMismatch mapsamples(+, x, rekey(x, 1, 1:4))
        @test_throws DimensionMismatch mapsamples(+, rekey(x, 2, 1:3), rekey(x, 2, 4:6))
        @test unwrap(mapsamples((a, _) -> a, x, z)) == data
    end
    @test_throws ArgumentError rekey(first(backends), 1, nothing)

    named = rebuild(last(backends), data, (Val(1), :channel => nothing))
    @test lookup(named, :channel) isa NoLookup
    @test axiskeys(rekey(named, 2, [:a, :b, :c]), 2) == [:a, :b, :c]

    x = first(backends)
    @test isnothing(axiskeys(@inferred(x[:, 2:3]), 2)) && isnothing(axiskeys(view(x, :, 2:3), 2))
    @test times(x[:, 2]) == t

    # Windows count samples without times; operations defined in time throw.
    windowed = [1.5 5.5 9.5; 3.5 7.5 11.5]
    @test tmean(data, 2; dim = 1) == windowed
    untimed = DimArray(data, (Ti(), Y(1:3)))
    r = tmean(untimed, 2)
    @test unwrap(r) == windowed && isnothing(axiskeys(r, 1))
    @test unwrap(mapsamples(+, last(backends), untimed)) == 2data
    @test_throws ArgumentError times(untimed)
    @test_throws ArgumentError tclip(untimed, 1, 2)
    @test_throws ArgumentError tderiv(untimed)
end
