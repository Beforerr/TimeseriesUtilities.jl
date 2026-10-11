@testitem "cadence" begin
    using Dates

    # https://github.com/JuliaSIMD/VectorizedStatistics.jl/issues/44#issue-3326201976
    t = Millisecond.(0:10000)
    @test cadence(t) == Millisecond(1)
    @test cadence(Millisecond.(0:20000)) == Millisecond(1)

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

    @test samplingrate(DateTime(2020) .+ Millisecond.(0:10:100)) == 100.0
end

@testitem "targmin, targmax, trange" begin
    using DimensionalData

    x = Timed([1.0 NaN NaN; NaN 5.0 NaN; 3.0 0.0 NaN], [10, 20, 30]; dim = 1)
    @test isequal(targmax(x), [30, 20, missing])
    @test isequal(targmin(x), [10, 30, missing])
    @test ismissing(targmax(Timed(Float64[], Int[])))
    @test isequal(tmaximum(x, 20), Timed([1.0 NaN NaN; 3.0 5.0 NaN], 0:20:20; dim = 1))

    t = Ti(1:10000)
    x = rand(t)
    @test trange(x) == (1, 10000)
    @test trange(x, x[Ti = 5:20], x[Ti = 10:30]) == (10, 20)
    @test isnothing(trange(x[Ti = 1:5], x[Ti = 6:9]))
end

@testitem "time grid" begin
    using Dates
    using Unitful
    using TimeseriesUtilities: _time_grid

    start_time = DateTime(2023, 1, 1, 0, 0, 0)
    end_time = DateTime(2023, 1, 1, 12, 0, 0)
    times = [start_time, start_time + Hour(3), start_time + Hour(6), end_time]
    @test _time_grid(times, Minute(30)) == start_time:Minute(30):end_time

    start_date = Date(2023, 1, 1)
    end_date = Date(2023, 1, 10)
    dates = [start_date, start_date + Day(3), start_date + Day(7), end_date]
    grid_daily = _time_grid(dates, Day(1))
    @test step(grid_daily) == Day(1)
    @test length(grid_daily) == 10

    @test _time_grid(times, 30u"minute") == _time_grid(times, Minute(30))
    @test _time_grid(times, Second(1)) == _time_grid(times, 1u"s") == _time_grid(times, 1u"Hz")

    @test _time_grid([start_time], Hour(1)) == [start_time]
end

@testitem "find_continuous_timeranges" begin
    using Dates
    using Chairmarks

    times = [
        Time(0, 0, 0), Time(1, 0, 0), Time(2, 0, 0),
        Time(6, 0, 0), Time(7, 0, 0), Time(8, 0, 0),
        Time(12, 0, 0), Time(13, 0, 0),
    ]

    @test find_continuous_timeranges(times, Hour(2)) ==
        [(Time(0, 0, 0), Time(2, 0, 0)), (Time(6, 0, 0), Time(8, 0, 0)), (Time(12, 0, 0), Time(13, 0, 0))]

    @test (@b find_continuous_timeranges($times, Hour(2))).allocs <= 2
end

@testitem "IntervalRange and tsplit" begin
    using Dates
    using Chairmarks

    # Test tsplit with integer n
    result = tsplit(1.0, 10, 3)
    @test collect(result) == [(1.0, 4.0), (4.0, 7.0), (7.0, 10.0)]

    result = tsplit(Date("2020-02-29"), Date("2022-03-01"), 4)
    @test result == [
        (Date("2020-02-29"), Date("2020-08-30")),
        (Date("2020-08-30"), Date("2021-03-01")),
        (Date("2021-03-01"), Date("2021-08-31")),
        (Date("2021-08-31"), Date("2022-03-01")),
    ]

    # Test tsplit with step size dt
    result = tsplit(0.0, 10.0, 2.0)
    @test length(result) == 5
    @test result[1] == (0.0, 2.0)
    @test result[5] == (8.0, 10.0)

    # Test tsplit with DateTime and Month period
    result = tsplit(DateTime("2021-01-01"), DateTime("2021-03-01"), Month)
    @test result == [
        (DateTime("2021-01-01"), DateTime("2021-02-01")),
        (DateTime("2021-02-01"), DateTime("2021-03-01")),
    ]

    eom = IntervalRange(Date("2021-01-31"), Date("2021-04-01"), Month(1))
    @test collect(eom) == [
        (Date("2021-01-31"), Date("2021-02-28")),
        (Date("2021-02-28"), Date("2021-03-31")),
        (Date("2021-03-31"), Date("2021-04-01")),
    ]
    @test eom[2] == (Date("2021-02-28"), Date("2021-03-31"))
    @test (@b collect(eom)).allocs <= 2

    quarter_split = IntervalRange(Date("2021-01-31"), Date("2021-08-01"), Quarter(1))
    @test collect(quarter_split) == [
        (Date("2021-01-31"), Date("2021-04-30")),
        (Date("2021-04-30"), Date("2021-07-31")),
        (Date("2021-07-31"), Date("2021-08-01")),
    ]

    year_split = IntervalRange(Date("2020-02-29"), Date("2022-03-01"), Year(1))
    @test year_split[3] == (Date("2022-02-28"), Date("2022-03-01"))
    @test_throws BoundsError year_split[4]

    # Test tsplit with DateTime and Day period
    result = tsplit(DateTime("2021-01-01"), DateTime("2021-01-05"), Day)
    @test result == [
        (DateTime("2021-01-01"), DateTime("2021-01-02")),
        (DateTime("2021-01-02"), DateTime("2021-01-03")),
        (DateTime("2021-01-03"), DateTime("2021-01-04")),
        (DateTime("2021-01-04"), DateTime("2021-01-05")),
    ]

    # Test tsplit with DateTime and Hour period
    result = tsplit(DateTime("2021-01-01T00:00:00"), DateTime("2021-01-01T03:00:00"), Hour)
    @test result == [
        (DateTime("2021-01-01T00:00:00"), DateTime("2021-01-01T01:00:00")),
        (DateTime("2021-01-01T01:00:00"), DateTime("2021-01-01T02:00:00")),
        (DateTime("2021-01-01T02:00:00"), DateTime("2021-01-01T03:00:00")),
    ]

    using JET
    @test_call tsplit(1.0, 10.0, 3)
    @test_call tsplit(0.0, 10.0, 2.0)
    @test_call tsplit(DateTime("2021-01-01"), DateTime("2021-03-01"), Month)
end

@testitem "TimeOffsets" begin
    using Dates
    using TimeseriesUtilities: TimeOffsets
    to = TimeOffsets(Day.(0:3), DateTime(2020, 1, 1))
    @test length(to) == 4
    @test to[[1, 4]] == DateTime.(2020, 1, [1, 4])
    to = TimeOffsets(1:10)
    @test length(to) == 10
    @test to[[1, 10]] == DateTime.(1970, 1, 1) .+ Second.([1, 10])
end
