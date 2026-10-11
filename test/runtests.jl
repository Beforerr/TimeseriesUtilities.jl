using TestItems, TestItemRunner

@run_package_tests

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(TimeseriesUtilities)
end

@testitem "api" begin
    using TimeseriesUtilities: axiskeys, rebuild
    @test_throws MethodError axiskeys([1, 2, 3], 1)
    @test_throws ArgumentError times(rand(2, 5))
    @test rebuild([1, 2, 3], 2:4, 1, 1:3) == 2:4
end
