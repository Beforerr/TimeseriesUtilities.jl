@testitem "cadence" begin
    using Dates
    using DimensionalData

    # https://github.com/JuliaSIMD/VectorizedStatistics.jl/issues/44#issue-3326201976
    t = Millisecond.(0:10000)
    @test cadence(t) == Millisecond(1)
    @test cadence(Millisecond.(0:20000)) == Millisecond(1)

    tdim = Ti(t)
    x = rand(tdim)
    @test cadence(x) == Millisecond(1)

    gapped = Millisecond.([0, 1, 2, 10, 11, 12])
    @test_logs (:warn, r"not approximately constant") cadence(gapped)
    @test cadence(gapped; check = false) == Millisecond(1)
    # Mixed cadence: gaps are not integer multiples → warn
    mixed = Millisecond.([0, 1, 2, 5, 6, 7])
    @test_logs (:warn, r"not approximately constant") cadence(mixed)

    # Robustness: 40% dropout (keeps 6/10 per cycle); diffs are [1,2,1,2,1,3,...]
    rng = 1:100
    all_t = Millisecond.(rng)
    kept = sort(all_t[filter(i -> mod(i * 7 + 3, 10) > 3, rng)])
    @test cadence(kept; check = false) == Millisecond(1)

    # Majority-gap dropout: pair+single pattern (period 7); 20 ones vs 39 threes
    # → median = 3, modal = 1
    kept2 = Millisecond.(vcat([[7k + 1, 7k + 2, 7k + 5] for k in 0:19]...))
    @test cadence(kept2; check = false) == Millisecond(1)

    # Large isolated burst gap (10 000× cadence) in otherwise regular data
    t_gap = vcat(Millisecond.(1:50), Millisecond.(10_001:10_050))
    @test cadence(t_gap; check = false) == Millisecond(1)

    # Realistic jitter: ±3% noise on 1 s cadence (Float64 timestamps)
    dt = 1.0
    t_jitter = dt .* (1:200) .+ 0.03dt .* (sin.(1:200))
    @test isapprox(cadence(t_jitter; check = false), dt; rtol = 0.01)

    # Mixed 2× cadence contamination: every 3rd point removed → 2Δt gaps
    base = collect(1:100)
    deleteat!(base, 3:3:99)
    @test cadence(base; check = false) == 1
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
    @test targmin(x) == t[argmin(x)]
    @test targmax(x) == t[argmax(x)]

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
    using Unitful

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

@testitem "find_continuous_timeranges" begin
    using AxisKeys
    using Dates
    using DimensionalData
    using Chairmarks

    # Create time series with gaps
    times = [
        Time(0, 0, 0), Time(1, 0, 0), Time(2, 0, 0), # Gap of 4 hours
        Time(6, 0, 0), Time(7, 0, 0), Time(8, 0, 0), # Gap of 4 hours
        Time(12, 0, 0), Time(13, 0, 0),
    ]

    # Find continuous ranges with max gap of 2 hours
    ranges = collect(ContinuousTimeRanges(times, Hour(2)))

    # Should find 3 continuous ranges
    @test length(ranges) == 3
    @test ranges[1] == (Time(0, 0, 0), Time(2, 0, 0))
    @test ranges[2] == (Time(6, 0, 0), Time(8, 0, 0))
    @test ranges[3] == (Time(12, 0, 0), Time(13, 0, 0))

    @test (@b ContinuousTimeRanges($times, Hour(2))).allocs == 0
    @test (@b find_continuous_timeranges($times, Hour(2))).allocs <= 2

    # Test with DimArray
    da = DimArray(rand(length(times)), (Ti(times),))
    ranges_da = find_continuous_timeranges(da, Hour(2))
    @test ranges_da == ranges

    ka = KeyedArray(rand(length(times)); time = times)
    ranges_ka = find_continuous_timeranges(ka, Hour(2))
    @test ranges_ka == ranges
end

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
