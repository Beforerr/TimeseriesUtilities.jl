using TestItems, TestItemRunner

@run_package_tests

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(TimeseriesUtilities)
end

@testitem "api" begin
    using TimeseriesUtilities: axiskeys, dims, rebuild_axes
    @test_throws MethodError axiskeys([1, 2, 3], 1)
    @test_throws MethodError dims([1, 2, 3], 1)
    @test rebuild_axes([1, 2, 3], 2:4, 1, 1:3) == 2:4
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
