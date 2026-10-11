@testitem "AxisKeys tgroupby" begin
    using AxisKeys
    using TimeseriesUtilities
    using Chairmarks

    ka = KeyedArray(10.0:10.0:50.0; time = 1.0:5.0)
    groups = tgroupby(ka, 2.0)

    @test length(groups) == 3
    @test groups[1] == [10.0]
    @test axiskeys(groups[2], 1) == [2.0, 3.0]
    @test axiskeys(groups[3], 1) == [4.0, 5.0]

    ka_large = KeyedArray(rand(1000); time = 1.0:1000.0)
    b = @b tgroupby($ka_large, 24.0)
    @info "tgroupby (1000 pts, window=24)" b
    @test b.allocs <= 270
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

@testitem "groupby_dynamic" begin
    using Dates
    using TimeseriesUtilities: GroupByDynamic, groupby_dynamic

    times = 1:1000
    dt = 24
    @test first.(collect(GroupByDynamic(times, dt))) == first(groupby_dynamic(times, dt))
    @test length(groupby_dynamic(times, dt)[1]) == 42
    times = Hour.(times)
    @test length(groupby_dynamic(times, Hour(dt))[1]) == 42
    times += DateTime(2010)
    @test length(groupby_dynamic(times, Hour(dt))[1]) == 42

    using JET
    @test_opt groupby_dynamic(times, Hour(dt))
    @test_call groupby_dynamic(times, Hour(dt))

    using Chairmarks
    verbose = true
    verbose && @info "groupby_dynamic" @b(groupby_dynamic($times, Hour($dt)))
    @test @b(groupby_dynamic($times, Hour($dt))).allocs <= 4
end
