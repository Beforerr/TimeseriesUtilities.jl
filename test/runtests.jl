using TestItems, TestItemRunner

@run_package_tests

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(TimeseriesUtilities)
end
