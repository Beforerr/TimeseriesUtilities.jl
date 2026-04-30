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
    using DataInterpolations: LinearInterpolation, ExtrapolationType
    using TimeseriesUtilities
    using TimeseriesUtilities: Tinterp

    times = [DateTime(2020, 1, 1), DateTime(2020, 1, 2), DateTime(2020, 1, 3)]
    da = DimArray(0:2, (Ti(times),))

    t = DateTime(2020, 1, 1, 12)
    @test tinterp(da, t) == tinterp(da, t; interp = Tinterp(LinearInterpolation))

    before = DateTime(2019, 12, 31)
    @test tinterp(da, before; extrapolation = true) ≈ -1.0
    @test tinterp(da, before; interp = Tinterp(LinearInterpolation), extrapolation = ExtrapolationType.Linear) ≈ -1.0
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
    @test_opt broken = true tsync(da1, da2, da3) # runtime dispatch
    @test_call tsync(da1, da2, da3)
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
