@testitem "twindows" begin
    using DimensionalData

    x = DimArray([1.0 10.0; 2.0 20.0; 3.0 30.0; 9.0 90.0], (Ti([0.0, 1.0, 2.5, 7.0]), Y(1:2)))
    ws, starts = twindows(x, 2.0; period = 3.0)
    @test starts === 0.0:2.0:6.0
    @test collect(TimeseriesUtilities.PointWindows(1:2, 0.5, 0.5)) == [(0.5, 1.5), (1.5, 2.5)]
    @test @inferred(TimeseriesUtilities.GridWindows(0.0f0:1.0f0:2.0f0, 1.0)[1]) === (0.0f0, 1.0)
    @test lookup(ws[2], Ti) == [2.5] && isempty(ws[3])
    # same windows as `tstat`
    @test reduce(vcat, [sum(parent(w); dims = 1) for w in ws]) == parent(tsum(x, 2.0; period = 3.0))
end

@testitem "window stops land on later starts" begin
    using Dates

    # from Jan 31, starts clamp to Feb 29, Mar 31, Apr 30
    x = Timed(ones(91), DateTime(2024, 1, 31) .+ Day.(0:90))
    @test parent(tsum(x, Month(1); period = Month(2), start_by = :datapoint))[2] == 61
    # `0.1 + 0.2` rounds above the start at 0.3
    y = Timed(ones(5), [0.1, 0.2, 0.3, 0.4, 0.5])
    @test parent(tsum(y, 0.1; period = 0.2, start_by = :datapoint)) == [2, 2, 2, 2, 1]
end

@testitem "float rounding does not move samples between windows" begin
    using Unitful

    for t in (-0.1 .+ (0:10) .* 0.1, Float32.(0:0.1:1))
        @test parent(tsum(Timed(ones(11), t), 0.1; start_by = :datapoint)) == ones(11)
    end
    # a sample at the rounding allowance of a stop on a later start
    t = [0.1, 0.2, prevfloat(0.3, 4), 0.4, 0.5]
    @test parent(tsum(Timed(ones(5), t), 0.1; period = 0.2, start_by = :datapoint)) == [2, 2, 2, 2, 1]
    # `mod` of a sample on the grid returns almost a whole step
    @test parent(tsum(Timed(ones(3), 0.7:0.1:0.9), 0.1)) == ones(3)
    # the last sample rounds to just below the end of the starts' range
    @test parent(tsum(Timed(ones(2), 0.3 .+ (0:1) .* 0.1), 0.1; start_by = :datapoint)) == ones(2)
    @test parent(tsum(Timed(ones(101), (0:0.1:10) .* u"s"), 0.1 * u"s")) == ones(101)
end

@testitem "rounding allowance leaves exact and infinite edges alone" begin
    # expanding window
    @test parent(tsmooth(sum, Timed(1.0:3.0, 0.0:2.0), (Inf, 0.5))) == [1.0, 3.0, 6.0]
    # integer epoch nanoseconds, beyond Float64 resolution
    @test first(twindows(Timed(ones(2), 1_700_000_000_000_000_123 .+ (0:1)), 1000)[2]) === 1_700_000_000_000_000_000
end

@testitem "tsmooth" begin
    using Dates
    using DimensionalData
    using TimeseriesUtilities

    data = [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0]
    times = [0.0, 1.0, 3.0, 6.0]

    # plain arrays smooth over sample indices
    @test tsmooth(data, 3; dim = 1) == [1.5 15.0; 2.0 20.0; 3.0 30.0; 3.5 35.0]
    @test tsmooth(data', 3; dim = 2) == [1.5 2.0 3.0 3.5; 15.0 20.0 30.0 35.0]
    @test parent(tsmooth(Timed(data, times; dim = 1), 2.0)) == [1.0 10.0; 1.5 15.0; 3.0 30.0; 4.0 40.0]
    @test parent(tsmooth(Timed(data', times; dim = 2), 2.0)) == [1.0 1.5 3.0 4.0; 10.0 15.0 30.0 40.0]
    # every backend smooths over its keys
    using AxisKeys
    @test tsmooth(KeyedArray(data; time = times, c = 1:2), 2.0) == [1.0 10.0; 1.5 15.0; 3.0 30.0; 4.0 40.0]

    @test parent(tsmooth(maximum, Timed(data, times; dim = 1), 2.0)) == [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0]

    explicit_left = parent(tsmooth(Timed(data, times; dim = 1), (2.0, 0.0)))
    @test isnan(explicit_left[1, 1])
    @test isnan(explicit_left[1, 2])

    t0 = DateTime(2020)
    datetimes = t0 .+ Second.([0, 1, 4, 5])
    da = DimArray([1.0, 2.0, 4.0, 8.0], (Ti(datetimes),))
    smoothed = tsmooth(da, Second(2))

    @test smoothed == [1.0, 1.5, 4.0, 6.0]
    @test parent(dims(smoothed, Ti)) == datetimes
    @test_throws ArgumentError tsmooth(da, 2)

    ms = t0 .+ Millisecond.(0:5)
    @test parent(tsmooth(Timed(1.0:6, ms), Millisecond(3))) == [1.5, 2.0, 3.0, 4.0, 5.0, 5.5]
    @test parent(tsmooth(Timed(1.0:6, ms), Millisecond(5))) == [2.0, 2.5, 3.0, 4.0, 4.5, 5.0]
    @test parent(tsmooth(Timed(1.0:6, Date(2020) .+ Day.(0:5)), Day(3))) == [1.5, 2.0, 3.0, 4.0, 5.0, 5.5]
end
