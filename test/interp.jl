@testitem "tinterp interpolation" begin
    using Dates, DimensionalData
    using TimeseriesUtilities

    # create a simple linearly increasing series
    times = [DateTime(2020, 1, 1), DateTime(2020, 1, 2), DateTime(2020, 1, 3)]
    da = DimArray(0:2, (Ti(times),))

    # single DateTime interpolation
    t1 = DateTime(2020, 1, 1, 12)
    res1 = tinterp(da, t1)
    @test res1 ≈ 0.5

    # multiple DateTime interpolation
    t2 = [DateTime(2020, 1, 1, 6), DateTime(2020, 1, 2, 18)]
    res2 = tinterp(da, t2)
    @test isa(res2, DimArray)
    @test res2 ≈ [0.25, 1.75]
    @test dims(res2, Ti).val == t2
    @test DimensionalData.span(dims(res2, Ti)) isa DimensionalData.Irregular
    t3 = DateTime(2020, 1, 1):Hour(6):DateTime(2020, 1, 2)
    @test DimensionalData.span(dims(tinterp(da, t3), Ti)) == DimensionalData.Regular(Hour(6))

    # create 3×2 series with numeric time dimension
    da3 = DimArray([1.0 4.0; 2.0 5.0; 3.0 6.0], (Ti(times), Y([10, 20])))

    # interpolate at two points
    res = tinterp(da3, t2)
    @test dims(res, Ti).val == t2
    @test res ≈ [1.25 4.25; 2.75 5.75]

    @test tinterp([0.0, 1.0, 2.0], [0.0, 1.0, 2.0], [0.5, 1.5]) ≈ [0.5, 1.5]
    @test tresample([0.0, 1.0, 2.0], [0.0, 1.0, 2.0], 1.0) ≈ [0.0, 1.0, 2.0]
end

@testitem "DataInterpolations compatibility" begin
    using Dates, DimensionalData
    using DataInterpolations: LinearInterpolation, ExtrapolationType, CubicSpline
    using TimeseriesUtilities
    using TimeseriesUtilities: Tinterp

    times = [DateTime(2020, 1, 1), DateTime(2020, 1, 2), DateTime(2020, 1, 3)]
    da = DimArray(0:2, (Ti(times),))

    t = DateTime(2020, 1, 1, 12)
    @test tinterp(da, t) == tinterp(da, t; interp = Tinterp(LinearInterpolation))

    before = DateTime(2019, 12, 31)
    @test tinterp(da, before; extrapolation = true) ≈ -1.0
    @test tinterp(da, before; interp = Tinterp(LinearInterpolation), extrapolation = ExtrapolationType.Linear) ≈ -1.0

    # CubicSpline on multi-dim array (per-channel fallback)
    t3 = [DateTime(2020, 1, i) for i in 1:5]
    da3 = DimArray([Float64(i) for i in 1:5, j in 1:3], (Ti(t3), Y(1:3)))
    t_new = [DateTime(2020, 1, 2, 12)]
    out = tinterp(da3, t_new; interp = Tinterp(CubicSpline))
    @test out isa DimArray
    @test size(out) == (1, 3)
    @test out[1, :] ≈ [2.5, 2.5, 2.5] atol = 0.1

end

@testitem "linear search hint" begin
    using TimeseriesUtilities: tinterp, LinearInterpolation
    # Repeated, backward and endpoint queries after a sorted run; nonlinear data so a
    # wrong segment changes the value. Reference: scalar path, which never uses the hint.
    t = cumsum(rand(50) .+ 0.1)
    u = sin.(t)
    q = [t[1]; sort(t[1] .+ (t[end] - t[1]) .* rand(200)); t[end]; t[3]; t[3]; t[20]]
    @test tinterp(u, t, q) ≈ LinearInterpolation(u, t).(q)
end

@testitem "FastInterpolations extension" begin
    using Dates, DimensionalData
    using FastInterpolations
    using TimeseriesUtilities

    t = collect(0.0:20.0)
    q = [0.5, 3.3, 19.9]
    A = randn(2, length(t), 3)  # time on dim 2, channels on dims 1 and 3
    for f in (cubic_interp, pchip_interp)  # `Series` path and per-channel path
        ref = [f(t, A[p, :, r], x) for p in 1:2, x in q, r in 1:3]
        @test tinterp(A, t, q; interp = f, dim = 2) ≈ ref
        @test tinterp(permutedims(A, (2, 1, 3)), t, q; interp = f, dim = 1) ≈ permutedims(ref, (2, 1, 3))
    end

    @test tinterp(t, t, [21.0]; interp = linear_interp, extrapolation = true) ≈ [21.0]

    times = DateTime(2020) .+ Hour.(0:3)
    da = DimArray([0.0, 1.0, 4.0, 9.0], (Ti(times),))
    @test tinterp(da, DateTime(2020, 1, 1, 1, 30); interp = cubic_interp) ≈ cubic_interp(0:3, parent(da), 1.5)
end

@testitem "AxisKeys tinterp" begin
    using AxisKeys
    using Dates
    using TimeseriesUtilities

    times = [DateTime(2020, 1, 1), DateTime(2020, 1, 2), DateTime(2020, 1, 3)]
    ka = KeyedArray(0.0:2.0; time = times)

    t1 = DateTime(2020, 1, 1, 12)
    @test tinterp(ka, t1) ≈ 0.5

    t2 = [DateTime(2020, 1, 1, 6), DateTime(2020, 1, 2, 18)]
    res = tinterp(ka, t2)
    @test res ≈ [0.25, 1.75]
    @test axiskeys(res, 1) == t2
    @test AxisKeys.dimnames(res, 1) == :time

    ka2 = KeyedArray(
        [1.0 4.0; 2.0 5.0; 3.0 6.0];
        time = times,
        component = [10, 20],
    )
    res2 = tinterp(ka2, t2)
    @test res2 ≈ [1.25 4.25; 2.75 5.75]
    @test axiskeys(res2, 1) == t2
    @test axiskeys(res2, 2) == [10, 20]
    @test AxisKeys.dimnames(res2, 2) == :component

    resampled = tresample(ka, Hour(12))
    @test axiskeys(resampled, 1) == DateTime(2020, 1, 1):Hour(12):DateTime(2020, 1, 3)
    @test resampled ≈ 0.0:0.5:2.0
end

@testitem "tsync" begin
    using Dates, DimensionalData
    include("./setup.jl")

    da1, da2, da3 = workload_interp_setup()
    a_sync, b_sync, c_sync = tsync(da1, da2, da3)

    # Check that all synchronized arrays have the same time dimension
    @test parent(dims(a_sync, Ti)) == parent(dims(b_sync, Ti)) == parent(dims(c_sync, Ti))

    # Check that the time range is the intersection of all input arrays
    @test dims(a_sync, Ti)[1] == DateTime(2020, 1, 2)
    @test dims(a_sync, Ti)[end] == DateTime(2020, 1, 3)

    # Check that values from the first and second array are preserved
    @test a_sync == [2, 3]
    @test b_sync == [10, 11]
    # The values should be interpolated at DateTime(2020, 1, 2) and DateTime(2020, 1, 3)
    expected_values = [
        5.5 9;
        6.5 11
    ]
    @test parent(c_sync) ≈ expected_values

    using JET
    @test_opt tsync(da1, da2, da3)
    @test_call tsync(da1, da2, da3)
end

@testitem "tfill_gaps" begin
    using Dates, DimensionalData
    using TimeseriesUtilities

    d(i, h = 0, m = 0) = DateTime(2020, 1, i, h, m)

    # jitter: count by rounding, fillers anchored to left original
    @test tfill_gaps([d(1), d(1, 23, 59), d(4, 0, 1)], Day(1)) ==
        [d(1), d(1, 23, 59), d(2, 23, 59), d(4, 0, 1)]

    # small gap filled, gap > max_gap left open
    @test tfill_gaps([d(1), d(3), d(10), d(11)], Day(1); max_gap = Day(2)) ==
        [d(1), d(2), d(3), d(10), d(11)]

    # (A, t): time on last dim, values preserved
    A = [1.0 2.0 4.0; 10.0 20.0 40.0]
    out, t_new = tfill_gaps(A, [0.0, 1.0, 3.0])
    @test t_new == [0.0, 1.0, 2.0, 3.0]
    @test isequal(out, [1.0 2.0 NaN 4.0; 10.0 20.0 NaN 40.0])

    # eltype: float precision kept, ints promoted to fit `fill`
    t = [0.0, 1.0, 3.0]
    @test eltype(first(tfill_gaps(Float32[1, 2, 4], t))) == Float32
    @test eltype(first(tfill_gaps([1, 2, 4], t))) == Float64
    @test isequal(first(tfill_gaps([1, 2, 4], t; fill = missing)), [1, 2, missing, 4])

    da = tfill_gaps(DimArray([1.0, 2.0, 4.0], (Ti([d(1), d(2), d(4)]),)))
    @test lookup(da, Ti) == [d(1), d(2), d(3), d(4)]
    @test isequal(parent(da), [1.0, 2.0, NaN, 4.0])
end

@testitem "tinterp_nans" begin
    using Dates, DimensionalData
    using TimeseriesUtilities

    # Create time series with NaN values in the middle
    times = [DateTime(2020, 1, 1), DateTime(2020, 1, 2), DateTime(2020, 1, 3), DateTime(2020, 1, 4), DateTime(2020, 1, 5)]
    data = [1.0, NaN, NaN, 4.0, 5.0]
    da = DimArray(data, (Ti(times),))

    # Interpolate NaN values
    result = tinterp_nans(da)

    @test result == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test result isa DimArray
    @test dims(result, Ti).val == times
    @test isnan(data[2])

    data2 = [1.0 10.0; NaN NaN; NaN 30.0; 4.0 40.0; 5.0 50.0]
    da2 = DimArray(data2, (Ti(times), Y([:a, :b])))
    result2 = tinterp_nans(da2)
    @test result2 ≈ [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0; 5.0 50.0]
    @test dims(result2, Ti).val == times
    @test dims(result2, Y).val == [:a, :b]
    @test isnan(data2[2, 1])
end
