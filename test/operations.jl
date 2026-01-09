@testitem "tsort, tclip, tview" begin
    using DimensionalData

    times = [3.0, 1.0, 2.0]
    t = Ti(times)
    y = Y(1:3)
    da = rand(t, y)
    result = tsort(da)
    @test parent(dims(result, Ti)) == sort(times)
    @test result[Ti(1)] == da[Ti(2)]
    @test result[Ti(3)] == da[Ti(1)]

    # tclip
    @test_throws "Cannot use an interval or `Between` with `Unordered`" tclip(da, 1.0, 2.0) == da[2:3, :]
    @test tclip(result, 1.0, 2.0) == parent(da[2:3, :])

    # tview
    @test_throws "Cannot use an interval or `Between` with `Unordered`" tview(da, 1.0, 2.0)
    @test tview(result, 1.0, 2.0) == parent(da[2:3, :])

    # Benchmark
    using Chairmarks
    verbose = false
    tclip_bench = @b tclip($result, 1.0, 2.0)
    tview_bench = @b tview($result, 1.0, 2.0)
    @test tclip_bench.allocs > tview_bench.allocs
    @test tclip_bench.time > tview_bench.time
    verbose && @info tclip_bench, tview_bench

    using JET
    @test_opt broken = true tsort(da) # TODO: `set` is not type-stable
    @test_call tsort(da)
    @test_opt tclip(result, 1.0, 2.0)
    @test_call tclip(result, 1.0, 2.0)
    @test_opt tview(result, 1.0, 2.0)
    @test_call tview(result, 1.0, 2.0)
end

@testitem "tselect" begin
    using DimensionalData

    # Test with numeric time dimension
    times = [1.0, 3.0, 5.0, 7.0, 9.0]

    @test tselect(times, 5.0) == 5.0
    @test tselect(times, 6.0) == 5.0
    @test tselect(times, 4.0) == 3.0

    values = [10.0, 20.0, 30.0, 40.0, 50.0]
    da = DimArray(values, (Ti(times),))

    # Test exact match
    @test tselect(da, 5.9) == tselect(da, 5.0, 0.5) == tselect(da, 5.2, 0.5) == tselect(da, 4.8, 0.5) == 30.0
    # Test closest match at edge of range
    @test tselect(da, 5.5, 0.5) == tselect(da, 4.5, 0.5) == 30.0

    # Test no match within range
    @test ismissing(tselect(da, 2.0, 0.5))
    @test ismissing(tselect(da, 8.0, 0.5))

    using JET
    @test_opt tselect(da, 5.0, 0.5)
    @test_call tselect(da, 5.0, 0.5)
end

@testitem "tsplit" begin
    using Dates

    # Test tsplit with integer n
    result = tsplit(1.0, 10.0, 3)
    @test length(result) == 3
    @test result[1] == (1.0, 4.0)
    @test result[2] == (4.0, 7.0)
    @test result[3] == (7.0, 10.0)

    # Test tsplit with step size dt
    result = tsplit(0.0, 10.0, 2.0)
    @test length(result) == 5
    @test result[1] == (0.0, 2.0)
    @test result[5] == (8.0, 10.0)

    # Test tsplit with DateTime and Month period
    result = tsplit(DateTime("2021-01-01"), DateTime("2021-03-01"), Month)
    @test length(result) == 2
    @test result[1] == (DateTime("2021-01-01"), DateTime("2021-02-01"))
    @test result[2] == (DateTime("2021-02-01"), DateTime("2021-03-01"))

    # Test tsplit with DateTime and Day period
    result = tsplit(DateTime("2021-01-01"), DateTime("2021-01-05"), Day)
    @test length(result) == 4
    @test result[1] == (DateTime("2021-01-01"), DateTime("2021-01-02"))
    @test result[2] == (DateTime("2021-01-02"), DateTime("2021-01-03"))
    @test result[3] == (DateTime("2021-01-03"), DateTime("2021-01-04"))
    @test result[4] == (DateTime("2021-01-04"), DateTime("2021-01-05"))

    # Test tsplit with DateTime and Hour period
    result = tsplit(DateTime("2021-01-01T00:00:00"), DateTime("2021-01-01T03:00:00"), Hour)
    @test length(result) == 3
    @test result[1] == (DateTime("2021-01-01T00:00:00"), DateTime("2021-01-01T01:00:00"))
    @test result[2] == (DateTime("2021-01-01T01:00:00"), DateTime("2021-01-01T02:00:00"))
    @test result[3] == (DateTime("2021-01-01T02:00:00"), DateTime("2021-01-01T03:00:00"))

    using JET
    @test_call tsplit(1.0, 10.0, 3)
    @test_call tsplit(0.0, 10.0, 2.0)
    @test_call tsplit(DateTime("2021-01-01"), DateTime("2021-03-01"), Month)
end

@testitem "tshift" begin
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    values = [10.0, 20.0, 30.0, 40.0, 50.0]

    @testset "DimensionalData" begin
        using DimensionalData
        da = DimArray(values, (Ti(times),))

        # Test shift by 2.0
        result = tshift(da, 2.0)
        @test parent(dims(result, Ti)) == [-1.0, 0.0, 1.0, 2.0, 3.0]

        # Test default shift with 2D array
        y = Y(1:3)
        da2d = rand(Ti(times), y)
        result2d = tshift(da2d)
        @test parent(dims(result2d, Ti)) == times .- 1.0
    end

    @testset "AxisKeys" begin
        using AxisKeys
        ka = KeyedArray(values; time = times)

        # Test shift by 2.0
        result = tshift(ka, 2.0)
        @test axiskeys(result, 1) == [-1.0, 0.0, 1.0, 2.0, 3.0]

        # Test with 2D array
        ka2d = KeyedArray(rand(5, 3); time = times, space = 1:3)
        result2d = tshift(ka2d)
        @test axiskeys(result2d, 1) == times .- 1.0
        @test axiskeys(result2d, 2) == 1:3
    end
end
