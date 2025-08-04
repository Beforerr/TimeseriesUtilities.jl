@testitem "tcross, tdot, tsproj, tnorm, tnorm_combine" begin
    using DimensionalData
    using LinearAlgebra
    t = Ti(1:10)
    y = Y(1:3)
    A = rand(t, y)
    B = rand(t, y)

    @test TimeseriesUtilities.proj([1, 2, 3], [4, 0, 0]) == [1, 0, 0]
    @test TimeseriesUtilities.proj!([1, 1, 1], [1, 2, 3], [4, 0, 0]) == [1, 0, 0]

    @test tcross(A, B)[1, :] == cross(A[1, :], B[1, :])
    @test tdot(A, B)[1] == dot(A[1, :], B[1, :])
    @test tsproj(A, B)[1] == sproj(A[1, :], B[1, :])
    @test tproj(A, B)[1, :] ≈ proj(A[1, :], B[1, :])
    @test toproj(A, B)[1, :] ≈ oproj(A[1, :], B[1, :])
    @test tnorm(A)[1] == norm(A[1, :])
    @test tnorm_combine(A)[1, :] ≈ [A[1, :]; norm(A[1, :])]

    using Chairmarks
    @test @b(tcross($A, $B)).allocs ≤ 2
    @test @b(tdot($A, $B)).allocs ≤ 2
    @test @b(tsproj($A, $B)).allocs ≤ 2
    @test @b(tproj($A, $B)).allocs ≤ 2
    @test @b(toproj($A, $B)).allocs ≤ 2
    @test @b(tnorm($A)).allocs ≤ 2

    using JET
    for f in (:tcross, :tdot, :tsproj, :tproj, :toproj)
        @eval @test_opt ignored_modules = (Base,) $f(A, B)
        @eval @test_call ignored_modules = (Base,) $f(A, B)
    end
end

@testitem "tderiv" begin
    using DimensionalData
    using Dates
    using TimeseriesUtilities.Unitful

    # Test with simple arrays and numeric times
    data = (1:5) .^ 2
    times = 1:5
    result = tderiv(data, times)
    expected = (1:4) .* 2 .+ 1  # derivative of x^2 is 2x +1
    @test result == expected

    # Test with 2D array along different dimensions
    data2d = [1.0 2.0; 4.0 8.0; 9.0 18.0; 16.0 32.0]  # [x^2 2x^2] for x = 1,2,3,4
    result2d = tderiv(data2d, 1:4; dim = 1)
    expected2d = [3.0 6.0; 5.0 10.0; 7.0 14.0]  # derivatives along time dimension
    @test result2d == expected2d

    # Test along second dimension
    data_dim2 = [1.0 4.0 9.0; 2.0 8.0 18.0]
    times_dim2 = [1.0, 2.0, 3.0]
    result_dim2 = tderiv(data_dim2, times_dim2; dim = 2)
    expected_dim2 = [3.0 5.0; 6.0 10.0]  # derivatives along second dimension
    @test result_dim2 == expected_dim2

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

    # Test with Date times (should return Unitful quantities)
    date_times = [Date(2020, 1, 1), Date(2020, 1, 2), Date(2020, 1, 3)]
    data_dates = [1.0, 4.0, 9.0]
    result_dates = tderiv(data_dates, date_times)
    @test result_dates[1] == 3.0u"d^-1"
    @test result_dates[2] == 5.0u"d^-1"

    # Test error handling
    @test_throws ArgumentError tderiv([1, 2, 3], [1, 2, 3]; dim = 2)  # dimension out of range

    # Test performance (should have minimal allocations)
    using Chairmarks
    A_perf = rand(Ti(1.0:100.0), Y(1:50))
    @test tderiv(A_perf) == tderiv(A_perf; lazy=true)
    # @info "benchmarks" b1 b2
    @test @b(tderiv($A_perf)).allocs ≤ 3
    @test @b(tderiv($A_perf; lazy=true)).allocs == 0
    @test @b(sum($tderiv($A_perf))).time > @b(sum($tderiv($A_perf; lazy=true))).time

    @info @b(tderiv($A_perf))

    # Test type inference
    using JET
    @test_opt ignored_modules = (Base,) tderiv(data, times)
    @test_call ignored_modules = (Base,) tderiv(data, times)
    @test_opt ignored_modules = (Base,) tderiv(A)
    @test_call ignored_modules = (Base,) tderiv(A)
end
