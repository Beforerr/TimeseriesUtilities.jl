@testitem "dropna" begin
    using DimensionalData

    t = Ti(2:5)
    y = Y([:x, :y])
    ds = DimStack(
        (
            a = DimArray([1.0, NaN, 3.0, 4.0], (t,)),
            b = DimArray([1.0 5.0; 2.0 NaN; 3.0 7.0; 4.0 8.0], (t, y)),
        )
    )

    a_clean = dropna(ds.a)
    @test a_clean == [1.0, 3.0, 4.0]
    @test a_clean.dims[1].val == [2, 4, 5]
    result = dropna(ds)
    @test dims(result, Ti).val == [2, 4, 5]
    @test result.a == [1.0, 3.0, 4.0]
    @test result.b == [1.0 5.0; 3.0 7.0; 4.0 8.0]
end

@testitem "tfill_gaps" begin
    using Dates, DimensionalData
    using TimeseriesUtilities: fill_times

    d(i, h = 0, m = 0) = DateTime(2020, 1, i, h, m)

    # jitter: count by rounding, fillers anchored to left original
    @test fill_times([d(1), d(1, 23, 59), d(4, 0, 1)], Day(1)) ==
        [d(1), d(1, 23, 59), d(2, 23, 59), d(4, 0, 1)]

    # small gap filled, gap > max_gap left open
    @test fill_times([d(1), d(3), d(10), d(11)], Day(1); max_gap = Day(2)) ==
        [d(1), d(2), d(3), d(10), d(11)]

    # values preserved
    A = [1.0 2.0 4.0; 10.0 20.0 40.0]
    out = tfill_gaps(Timed(A, [0.0, 1.0, 3.0]; dim = 2))
    @test times(out) == [0.0, 1.0, 2.0, 3.0]
    @test isequal(parent(out), [1.0 2.0 NaN 4.0; 10.0 20.0 NaN 40.0])

    # eltype: float precision kept, ints promoted to fit `fill`
    t = [0.0, 1.0, 3.0]
    @test eltype(parent(tfill_gaps(Timed(Float32[1, 2, 4], t)))) == Float32
    @test eltype(parent(tfill_gaps(Timed([1, 2, 4], t)))) == Float64
    @test isequal(parent(tfill_gaps(Timed([1, 2, 4], t); fill = missing)), [1, 2, missing, 4])

    da = tfill_gaps(DimArray([1.0, 2.0, 4.0], (Ti([d(1), d(2), d(4)]),)))
    @test lookup(da, Ti) == [d(1), d(2), d(3), d(4)]
    @test isequal(parent(da), [1.0, 2.0, NaN, 4.0])
end

@testitem "tinterp_na" begin
    using Dates, DimensionalData
    using TimeseriesUtilities

    # Create time series with NaN values in the middle
    times = [DateTime(2020, 1, 1), DateTime(2020, 1, 2), DateTime(2020, 1, 3), DateTime(2020, 1, 4), DateTime(2020, 1, 5)]
    data = [1.0, NaN, NaN, 4.0, 5.0]
    da = DimArray(data, (Ti(times),))

    # Interpolate NaN values
    result = tinterp_na(da)

    @test result == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test result isa DimArray
    @test dims(result, Ti).val == times
    @test isnan(data[2])

    data2 = [1.0 10.0; NaN NaN; NaN 30.0; 4.0 40.0; 5.0 50.0]
    da2 = DimArray(data2, (Ti(times), Y([:a, :b])))
    result2 = tinterp_na(da2)
    @test result2 ≈ [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0; 5.0 50.0]
    @test dims(result2, Ti).val == times
    @test dims(result2, Y).val == [:a, :b]
    @test isnan(data2[2, 1])
end
