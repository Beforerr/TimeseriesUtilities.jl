@testitem "proj" begin
    @test TimeseriesUtilities.proj([1, 2, 3], [4, 0, 0]) == [1, 0, 0]
end

@testitem "tderiv" begin
    using DimensionalData
    using Dates
    using Unitful

    # Test with simple arrays and numeric times
    data = Timed((1:5) .^ 2, 1:5)
    result = tderiv(data)
    @test parent(result) == (1:4) .* 2 .+ 1  # derivative of x^2 is 2x +1
    @test times(result) == 1:4

    data2d = [1.0 2.0; 4.0 8.0; 9.0 18.0; 16.0 32.0]  # [x^2 2x^2] for x = 1,2,3,4
    @test parent(tderiv(Timed(data2d, 1:4; dim = 1))) == [3.0 6.0; 5.0 10.0; 7.0 14.0]
    @test parent(tderiv(Timed(permutedims(data2d), 1.0:4; dim = 2))) == [3.0 5.0 7.0; 6.0 10.0 14.0]
    @test parent(tderiv(Timed((@view data2d[2:4, :]), 2:4; dim = 1))) == [5.0 10.0; 7.0 14.0]

    # Test with DimensionalData arrays
    t = Ti(1.0:1.0:5.0)
    y = Y(1:3)
    A = DimArray([i^2 for i in 1:5, j in 1:3], (t, y))

    # Test default behavior (should use time dimension)
    result_dd = tderiv(A)
    @test tderiv(A; dim = 1) == result_dd
    @test size(result_dd) == (4, 3)  # one less in time dimension
    @test result_dd[1, 1] == 3.0  # derivative at first point
    @test result_dd[2, 1] == 5.0  # derivative at second point

    # Test with Date times (returns values per second)
    date_times = [Date(2020, 1, 1), Date(2020, 1, 2), Date(2020, 1, 3)]
    data_dates = [1.0, 4.0, 9.0]
    @test parent(tderiv(Timed(data_dates, date_times))) == [3.0, 5.0] ./ 86400

    using Chairmarks
    A_perf = rand(Ti(1.0:100.0), Y(1:50))
    @test tderiv(A_perf) == tderiv(A_perf; lazy = true)
    @test @b(tderiv($A_perf)).allocs ≤ 3
    @test @b(tderiv($A_perf; lazy = true)).allocs == 0 broken = VERSION < v"1.11"

    # Test type inference
    using JET
    @test_opt ignored_modules = (Base,) tderiv(data)
    @test_call ignored_modules = (Base,) tderiv(data)
    @test_opt tderiv(A)
    @test_call ignored_modules = (Base,) tderiv(A)
end

@testitem "tsubtract" begin
    x = Timed([1.0 10.0 100.0; 3.0 30.0 300.0], [0.0, 1.0, 2.0]; dim = 2)
    @test parent(tsubtract(first, x)) == [0.0 9.0 99.0; 0.0 27.0 297.0]
    @test parent(tsubtract(sum, Timed([1.0, 2.0, 3.0], 1:3))) == [-5.0, -4.0, -3.0]
    # running median over 3 samples leaves the spike
    using Statistics: median
    @test parent(tsubtract(median, Timed([1.0, 2.0, 10.0, 4.0, 5.0], 1:5), 3)) == [-0.5, 0.0, 6.0, -1.0, 0.5]
end
