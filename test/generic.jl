@testitem "Generic: tsort, tclip, tview" begin
    times = [3.0, 1.0, 2.0, 4.0, 5.0]
    data = [30.0, 10.0, 20.0, 40.0, 50.0]

    # tsort
    sorted_data, sorted_times = tsort(data, times)
    @test sorted_times == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test sorted_data == [10.0, 20.0, 30.0, 40.0, 50.0]

    # tclip
    clipped_data, clipped_times = tclip(sorted_data, sorted_times, 2.0, 4.0)
    @test clipped_times == [2.0, 3.0, 4.0]
    @test clipped_data == [20.0, 30.0, 40.0]

    # tview
    viewed_data, viewed_times = tview(sorted_data, sorted_times, 2.0, 4.0)
    @test viewed_times == [2.0, 3.0, 4.0]
    @test viewed_data == [20.0, 30.0, 40.0]

    # 2D data
    data2d = [30.0 300.0; 10.0 100.0; 20.0 200.0; 40.0 400.0; 50.0 500.0]
    sorted2d, sorted_t = tsort(data2d, times)
    @test sorted_t == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test sorted2d == [10.0 100.0; 20.0 200.0; 30.0 300.0; 40.0 400.0; 50.0 500.0]
end

@testitem "Generic: tselect" begin
    times = [1.0, 3.0, 5.0, 7.0, 9.0]
    values = [10.0, 20.0, 30.0, 40.0, 50.0]

    @test tselect(times, 5.0) == 5.0
    @test tselect(times, 6.0) == 5.0

    @test tselect(values, times, 5.0) == 30.0
    @test tselect(values, times, 6.0) == 30.0

    @test tselect(values, times, 5.0, 0.5) == 30.0
    @test ismissing(tselect(values, times, 2.0, 0.5))
end

@testitem "Generic: tmask!" begin
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    data = [10.0, 20.0, 30.0, 40.0, 50.0]

    masked = tmask(data, times, 2.0, 4.0)
    @test masked[1] == 10.0
    @test isnan(masked[2])
    @test isnan(masked[3])
    @test isnan(masked[4])
    @test masked[5] == 50.0
end

@testitem "Generic: tshift" begin
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    values = [10.0, 20.0, 30.0, 40.0, 50.0]

    _, shifted_times = tshift(values, times, 2.0)
    @test shifted_times == [-1.0, 0.0, 1.0, 2.0, 3.0]

    _, default_shifted = tshift(values, times)
    @test default_shifted == [0.0, 1.0, 2.0, 3.0, 4.0]
end

@testitem "Generic: tderiv" begin
    data = Float64.((1:5) .^ 2)
    times = Float64.(1:5)
    result = tderiv(data, times)
    expected = Float64.((1:4) .* 2 .+ 1)
    @test result == expected

    # 2D
    data2d = [1.0 2.0; 4.0 8.0; 9.0 18.0; 16.0 32.0]
    result2d = tderiv(data2d, Float64.(1:4); dim=1)
    expected2d = [3.0 6.0; 5.0 10.0; 7.0 14.0]
    @test result2d == expected2d
end

@testitem "Generic: tstat / tmean / tmedian" begin
    using Statistics

    data1d = [1.0, 2.0, 3.0, 4.0]
    @test tstat(sum, data1d) == 10.0
    @test tmean(data1d) ≈ mean(data1d)
    @test tmedian(data1d) ≈ median(data1d)

    data2d = [1.0 5.0; 2.0 6.0; 3.0 7.0; 4.0 8.0]
    @test tstat(sum, data2d; dim=1) == [10.0, 26.0]
end

@testitem "Generic: tnorm, tdot, tcross" begin
    using LinearAlgebra

    A = rand(10, 3)
    B = rand(10, 3)

    @test tnorm(A; dim=1)[1] == norm(A[1, :])
    @test tdot(A, B; dim=1)[1] == dot(A[1, :], B[1, :])
    @test tcross(A, B; dim=1)[1, :] == cross(A[1, :], B[1, :])
end

@testitem "Generic: tinterp" begin
    using DataInterpolations

    old_times = [1.0, 2.0, 3.0, 4.0, 5.0]
    data = [10.0, 20.0, 30.0, 40.0, 50.0]
    new_times = [1.5, 2.5, 3.5, 4.5]

    result = tinterp(data, old_times, new_times; dim=1)
    @test result ≈ [15.0, 25.0, 35.0, 45.0]

    # 2D
    data2d = [10.0 100.0; 20.0 200.0; 30.0 300.0]
    result2d = tinterp(data2d, [1.0, 2.0, 3.0], [1.5, 2.5]; dim=1)
    @test result2d ≈ [15.0 150.0; 25.0 250.0]
end

@testitem "Generic: tgroupby" begin
    times = 1:100
    data = Float64.(1:100)

    groups = tgroupby(data, times, 25)
    @test length(groups) == 4
    @test groups[1] == Float64.(1:25)
end

@testitem "Generic: smooth" begin
    data = ones(20)
    result = smooth(data, 5; dim=1)
    @test all(result .≈ 1.0)
end

@testitem "Generic: tsubtract" begin
    data = [1.0 2.0; 3.0 4.0; 5.0 6.0]
    result = tsubtract(data; dim=1)
    # median along dim=1 is [3.0, 4.0], so result = data .- [3.0 4.0]
    @test result ≈ [-2.0 -2.0; 0.0 0.0; 2.0 2.0]
end

@testitem "Generic: tnorm_combine" begin
    using LinearAlgebra
    A = [1.0 0.0; 0.0 1.0; 1.0 1.0]
    result = tnorm_combine(A; dim=1)
    @test size(result) == (3, 3)
    @test result[1, 3] ≈ norm([1.0, 0.0])
    @test result[3, 3] ≈ norm([1.0, 1.0])
end

@testitem "Generic: tinterp_nans" begin
    data = [1.0, NaN, NaN, 4.0, 5.0]
    times = Float64.(1:5)
    result = tinterp_nans(data, times; dim=1)
    @test result[1] == 1.0
    @test result[4] == 4.0
    @test result[5] == 5.0
    @test result[2] ≈ 2.0
    @test result[3] ≈ 3.0
end
