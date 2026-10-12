@testitem "twindows" begin
    using DimensionalData

    x = DimArray([1.0 10.0; 2.0 20.0; 3.0 30.0; 9.0 90.0], (Ti([0.0, 1.0, 2.5, 7.0]), Y(1:2)))
    ws, starts = twindows(x, 2.0; period = 3.0)
    @test starts === 0.0:2.0:6.0
    @test collect(TimeseriesUtilities.PointWindows(1:2, 0.5, 0.5)) == [(0.5, 1.5), (1.5, 2.5)]
    @test lookup(ws[2], Ti) == [2.5] && isempty(ws[3])
    # same windows as `tstat`
    @test reduce(vcat, [sum(parent(w); dims = 1) for w in ws]) == parent(tsum(x, 2.0; period = 3.0))
end

@testitem "smooth" begin
    using Dates
    using DimensionalData
    using TimeseriesUtilities

    data = [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0]
    times = [0.0, 1.0, 3.0, 6.0]

    # plain arrays smooth over sample indices
    @test smooth(data, 3; dim = 1) == [1.5 15.0; 2.0 20.0; 3.0 30.0; 3.5 35.0]
    @test smooth(data', 3; dim = 2) == [1.5 2.0 3.0 3.5; 15.0 20.0 30.0 35.0]
    @test parent(smooth(Timed(data, times; dim = 1), 2.0)) == [1.0 10.0; 1.5 15.0; 3.0 30.0; 4.0 40.0]
    @test parent(smooth(Timed(data', times; dim = 2), 2.0)) == [1.0 1.5 3.0 4.0; 10.0 15.0 30.0 40.0]
    # every backend smooths over its keys
    using AxisKeys
    @test smooth(KeyedArray(data; time = times, c = 1:2), 2.0) == [1.0 10.0; 1.5 15.0; 3.0 30.0; 4.0 40.0]

    @test parent(smooth(maximum, Timed(data, times; dim = 1), 2.0)) == [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0]

    explicit_left = parent(smooth(Timed(data, times; dim = 1), (2.0, 0.0)))
    @test isnan(explicit_left[1, 1])
    @test isnan(explicit_left[1, 2])

    t0 = DateTime(2020)
    datetimes = t0 .+ Second.([0, 1, 4, 5])
    da = DimArray([1.0, 2.0, 4.0, 8.0], (Ti(datetimes),))
    smoothed = smooth(da, Second(2))

    @test smoothed == [1.0, 1.5, 4.0, 6.0]
    @test parent(dims(smoothed, Ti)) == datetimes
    @test_throws ArgumentError smooth(da, 2)

    ms = t0 .+ Millisecond.(0:5)
    @test parent(smooth(Timed(1.0:6, ms), Millisecond(3))) == [1.5, 2.0, 3.0, 4.0, 5.0, 5.5]
    @test parent(smooth(Timed(1.0:6, ms), Millisecond(5))) == [2.0, 2.5, 3.0, 4.0, 4.5, 5.0]
    @test parent(smooth(Timed(1.0:6, Date(2020) .+ Day.(0:5)), Day(3))) == [1.5, 2.0, 3.0, 4.0, 5.0, 5.5]
end
