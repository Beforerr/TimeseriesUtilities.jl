# References
# - https://www.mathworks.com/help/matlab/ref/isoutlier.html
# - https://www.mathworks.com/help/matlab/ref/filloutliers.html

@testitem "outlier_mask, despike" begin
    using LinearAlgebra
    A = [57, 59, 60, 100, 59, 58, 57, 58, 300, 61, 62, 60, 62, 58, 57]
    @test findall(outlier_mask(A)) == [4, 9]
    @test findall(isnan, despike(A)) == [4, 9]

    x = -2π:0.1:2π
    A = sin.(x)
    A[47] = 0
    A[60] = NaN
    @test findall(outlier_mask(A, 5)) == [47]
    # window is in time units, not samples
    @test parent(outlier_mask(Timed(A, 10 .* eachindex(A)), 50)) == outlier_mask(A, 5)

    A = rand(5, 5) + 200I
    @test sum(diag(outlier_mask(A; dim = 2))) == 5

    @test_throws ArgumentError despike!([1, 2, 100, 2, 1])
end

@testitem "rolling outlier_mask matches per-window statistics" begin
    using Random, Statistics
    function brute(x, t, (before, after); k = 3)
        return map(eachindex(x)) do i
            w = filter(!isnan, x[t[i] - before .<= t .< t[i] + after])
            isempty(w) && return false
            med = median(w)
            abs(x[i] - med) > k * 1.4826 * median(abs.(w .- med))
        end
    end
    Random.seed!(1)
    for _ in 1:500
        n = rand(1:60)
        x = round.(randn(n); digits = 1) # ties
        x[rand(1:n, rand(0:(n ÷ 3)))] .= NaN
        x[rand(1:n, 3)] .= 50
        t = cumsum(rand(1:8, n)) # uneven spacing, empty windows, samples on window edges
        window = (rand(0:12), rand(0:12))
        @test parent(outlier_mask(Timed(x, t), window)) == brute(x, t, window)
    end
end
