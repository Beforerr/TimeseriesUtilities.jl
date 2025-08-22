@testitem "resolution" begin
    using Dates
    using DimensionalData

    # https://github.com/JuliaSIMD/VectorizedStatistics.jl/issues/44#issue-3326201976
    t = Millisecond.(0:10000)
    @test resolution(t) == Millisecond(1)
    @test resolution(Millisecond.(0:20000)) == Millisecond(1)

    tdim = Ti(t)
    x = rand(tdim)
    @test resolution(x) == Millisecond(1)

    using Chairmarks
    using TimeseriesUtilities: median, nanmedian, vmedian
    t = Millisecond.(0:3:100000)
    dt = diff(t)
    
    @info "Benchmarking median implementations:"
    f1(dt) = median(reinterpret(Int, dt))
    f2(dt) = vmedian(reinterpret(Int, dt))
    f3(dt) = nanmedian(reinterpret(Int, dt))
    b1, b2, b3 = @b f1(dt), f2(dt), f3(dt)
    @test f1(dt) == f2(dt) == f3(dt)
    # Find fastest implementation
    times = [b1.time, b2.time, b3.time]
    names = ["Statistics.median", "VectorizedStatistics.vmedian", "NaNStatistics.nanmedian"]
    fastest_idx = argmin(times)
    @info "Times: $(names .=> times)"
    @info "Fastest median implementation: $(names[fastest_idx])"
end


@testitem "tmin, tmax, timerange" begin
    using Chairmarks
    using Dates
    using DimensionalData

    t = Ti(1:10000)
    x = rand(t)
    @test tminimum(x) == minimum(t)
    @test tmaximum(x) == maximum(t)
    @test timerange(x) == extrema(t)

    verbose = false

    for T in (Int, Date, DateTime)
        ts = T.(collect(1:10000))
        @test timerange(ts) == extrema(ts)
        b1 = @b timerange($ts)
        b2 = @b extrema($ts)
        if b1.time < b2.time
            verbose && @info "Acceleration ratio: $(b2.time / b1.time)"
        else
            @info "Deceleration ratio: $(b1.time / b2.time)"
        end
    end
end

@testitem "time_grid" begin
    using Dates
    using DimensionalData
    using TimeseriesUtilities.Unitful

    # Test with DateTime data
    start_time = DateTime(2023, 1, 1, 0, 0, 0)
    end_time = DateTime(2023, 1, 1, 12, 0, 0)
    times = [start_time, start_time + Hour(3), start_time + Hour(6), end_time]
    grid_30min = time_grid(times, Minute(30))
    @test first(grid_30min) == start_time
    @test last(grid_30min) == end_time
    @test step(grid_30min) == Minute(30)
    @test length(grid_30min) == 25  # Every 30 minutes for 12 hours

    # Test with Date data
    start_date = Date(2023, 1, 1)
    end_date = Date(2023, 1, 10)
    dates = [start_date, start_date + Day(3), start_date + Day(7), end_date]
    grid_daily = time_grid(dates, Day(1))
    @test step(grid_daily) == Day(1)
    @test length(grid_daily) == 10  # 10 days total

    # Test with Unitful
    @test_throws MethodError start_time:1u"hr":end_time
    @test grid_30min == time_grid(times, 30u"minute")
    @test time_grid(times, Second(1)) == time_grid(times, 1u"s") == time_grid(times, 1u"Hz")

    # Test with DimensionalData array
    t = Ti(times)
    x = rand(t)
    @test time_grid(x, Hour(2)) == time_grid(times, Hour(2))

    # Test edge case: single time point
    single_time = [DateTime(2023, 1, 1)]
    grid_single = time_grid(single_time, Hour(1))
    @test length(grid_single) == 1
    @test first(grid_single) == last(grid_single)
end
