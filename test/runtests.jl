using TestItems, TestItemRunner

@run_package_tests

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(TimeseriesUtilities)
end

@testitem "README covers exports" begin
    readme = read(joinpath(pkgdir(TimeseriesUtilities), "README.md"), String)
    @test filter(n -> !occursin(Regex("`\\Q$n\\E(?![\\w!])"), readme), setdiff(names(TimeseriesUtilities), [:TimeseriesUtilities])) == Symbol[]
end

@testitem "api" begin
    using TimeseriesUtilities: axiskeys, rebuild, resolve_dim, tdimnum, unwrap
    using DimensionalData: Ti
    @test isnothing(axiskeys([1, 2, 3], 1))
    @test_throws ArgumentError times([1, 2, 3])
    @test isnothing(tdimnum([1, 2])) && resolve_dim([1, 2], nothing) == 1
    @test_throws ArgumentError times(rand(2, 5))
    @test rebuild([1, 2, 3], 2:4, 1, 1:3) == 2:4
    # a container without `rebuild` must not silently drop its keys
    struct Wrapped{T} <: AbstractVector{T}
        data::Vector{T}
    end
    TimeseriesUtilities.unwrap(w::Wrapped) = w.data
    @test_throws ArgumentError rebuild(Wrapped([1, 2]), [3, 4])
    # a keyed vector has a time dimension only if it declares one
    @test_throws ArgumentError resolve_dim(Wrapped([1, 2]), nothing)
    @test unwrap(Ti(view([1, 2, 3, 4], 2:3))) == [2, 3]
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
    @test hash(x) == hash(Timed(copy(A), copy(t); dim = 2))
    @test isequal(Timed([NaN], [1]), Timed([NaN], [1]))
    data, times_ = x
    @test data === A && times_ === t
    @test_throws MethodError collect(x)
    @test TimeseriesUtilities.tdimnum(TimeseriesUtilities.rebuild(x, A', (Val(2), Val(1)))) == 1
end

@testitem "Timed keys varying in time" begin
    using TimeseriesUtilities: axiskeys, rebuild, rekey, unwrap
    data = reshape(1.0:12.0, 4, 3)
    t = 10:10:40
    E = [10.0 20 30; 11 21 31; 12 22 32; 13 23 33]  # (time, energy)
    x = Timed{1}(data, (t, E))
    xt = Timed{2}(permutedims(data), (permutedims(E), t))
    @test_throws DimensionMismatch Timed{1}(data, (t, permutedims(E)))

    @test axiskeys(x[2:3, 1:2], 2) == E[2:3, 1:2]
    @test axiskeys(xt[1:2, 2:3], 1) == permutedims(E)[1:2, 2:3]
    @test axiskeys(tshift(x), 2) === E
    @test isnothing(axiskeys(rebuild(x, data, 1, 1:4), 2))
    @test axiskeys(rebuild(x, permutedims(data), (Val(2), Val(1))), 1) == permutedims(E)

    @test unwrap(mapsamples(+, x, xt)) == 2data
    @test unwrap(mapsamples(+, rekey(x, 2, repeat(E[1:1, :], 4)), rekey(xt, 1, E[1, :]))) == 2data
    @test_throws DimensionMismatch mapsamples(+, x, rekey(xt, 1, E[1, :]))
end
