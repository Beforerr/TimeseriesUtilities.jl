using TestItems, TestItemRunner

@run_package_tests

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(TimeseriesUtilities)
end

@testitem "api" begin
    using TimeseriesUtilities: axiskeys, rebuild, unwrap
    @test_throws ArgumentError axiskeys([1, 2, 3], 1)
    @test_throws ArgumentError times(rand(2, 5))
    @test rebuild([1, 2, 3], 2:4, 1, 1:3) == 2:4
    # a container without `rebuild` must not silently drop its keys
    struct Wrapped{T} <: AbstractVector{T}
        data::Vector{T}
    end
    TimeseriesUtilities.unwrap(w::Wrapped) = w.data
    @test_throws ArgumentError rebuild(Wrapped([1, 2]), [3, 4])
end

@testitem "Timed" begin
    A = [1.0 2.0 3.0; 4.0 5.0 6.0]
    t = [10, 20, 30]
    x = Timed(A, t; dim = 2)
    @test_throws DimensionMismatch Timed(A, t[1:2]; dim = 2)
    @test_throws ArgumentError Timed(A, t)
    # indexing keeps the keys of kept axes and tracks where time went
    @test x[1, :] isa Timed{1} && times(x[1, :]) == t
    @test x[:, 2] == [2.0, 5.0]
    @test times(view(x, :, 2:3)) == [20, 30] && times(view(x, :, 2:3)) isa SubArray
    @test x[end, :] == Timed(A[2, :], t)
    @test Timed([1 2; 3 4], 1:2; dim = 1) != Timed([1 2; 3 4], 1:2; dim = 2)
    @test x == Timed(A, t; dim = Int32(2))
    data, times_ = x
    @test data === A && times_ === t
    @test_throws MethodError collect(x)
    @test TimeseriesUtilities.tdimnum(TimeseriesUtilities.rebuild(x, A', (Val(2), Val(1)), (t, 1:2))) == 1
end
