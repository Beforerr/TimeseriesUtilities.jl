@testitem "tcross, tdot, tsproj, tnorm" begin
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

    using Chairmarks
    @test @b(tcross($A, $B)).allocs ≤ 2
    @test @b(tdot($A, $B)).allocs ≤ 2
    @test @b(tsproj($A, $B)).allocs ≤ 2
    @test @b(tproj($A, $B)).allocs ≤ 2
    @test @b(toproj($A, $B)).allocs ≤ 2
    @test @b(tnorm($A)).allocs ≤ 2

    using JET
    for f in (:tcross, :tdot, :tsproj, :tproj, :toproj)
        @eval @test_opt ignored_modules=(Base,) $f(A, B)
        @eval @test_call ignored_modules=(Base,) $f(A, B)
    end
end
