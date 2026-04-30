@testitem "groupby_dynamic" begin
    using Dates
    using TimeseriesUtilities: groupby_dynamic

    times = 1:1000
    dt = 24
    @test length(groupby_dynamic(times, dt)[1]) == 42
    times = Hour.(times)
    @test length(groupby_dynamic(times, Hour(dt))[1]) == 42
    times += DateTime(2010)
    @test length(groupby_dynamic(times, Hour(dt))[1]) == 42

    using JET
    @test_opt groupby_dynamic(times, Hour(dt))
    @test_call groupby_dynamic(times, Hour(dt))

    using Chairmarks
    verbose = true
    verbose && @info "groupby_dynamic" @b(groupby_dynamic($times, Hour($dt)))
end

@testitem "timeseries statistics" begin
    using Dates
    using Statistics
    using Chairmarks
    using TimeseriesUtilities.NaNStatistics

    @testset "tmean" begin
        t0 = DateTime(2011)
        offsets = Millisecond.(0:9)
        ts = t0 .+ offsets
        @test_throws InexactError mean(offsets)
        @test_throws MethodError mean(ts)
        @test tmean(ts) == t0 + Microsecond(4500)
        @test_throws MethodError tmean([1.0, 2.0, 3.0], 1.0)
    end

    @testset "DimensionalData" begin
        using DimensionalData


        t = Ti(Millisecond.(0:3))
        y = Y(1:2)
        da1 = rand(t)
        da2 = rand(t, y)

        @test tmean(da1) == mean(da1)
        @test tmean(da1, Millisecond(2)) == [mean(da1[1:2]), mean(da1[3:4])]
        @test tmean(da2) == vec(mean(da2, dims = 1))
        @test tmean(da2, Millisecond(2)) == [mean(parent(da2)[1:2, :], dims = 1); mean(parent(da2)[3:4, :], dims = 1)]

        # tmedian
        @test tmedian(da1) == median(da1)
        @test tmedian(da1, Millisecond(2)) == [median(da1[1:2]), median(da1[3:4])]
        @test tmedian(da2) == vec(median(da2, dims = 1))

        # tsum, tvar, tstd, tsem
        @test tsum(da2) == vec(sum(da2, dims = 1))
        @test tvar(da1) == var(da1)
        @test tstd(da1) == std(da1)
        @test tsem(da1) == nansem(da1)

        # DimStack
        @test tmean(DimStack((da1, da2))) == (; layer1 = tmean(da1), layer2 = tmean(da2))

        @testset "output coordinate" begin
            grouped2 = tmean(da2, Millisecond(2))
            @test grouped2.dims[1].val == Millisecond.([0, 2])
            @test grouped2.dims[2].val == 1:2
        end

        verbose = false
        da_bench1 = rand(Ti(1:1000))
        verbose && @info "tmean" @b(tmean($da_bench1))
        da_bench = DimArray(rand(1000, 3), (Ti(1:1000), Y(1:3)))
        verbose && @info "tmean" @b(tmean($da_bench, 10))
    end

    @testset "AxisKeys" begin
        using AxisKeys

        times = Millisecond.(0:3)
        ka1 = KeyedArray([1.0, 2.0, 3.0, 4.0]; time = times)
        ka2 = KeyedArray(
            [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0];
            time = times,
            component = [1, 2],
        )

        grouped1 = tmean(ka1, Millisecond(2))
        @test grouped1 == [1.5, 3.5]
        @test AxisKeys.dimnames(grouped1, 1) == :time

        grouped2 = tmean(ka2, Millisecond(2))
        @test grouped2 == [1.5 15.0; 3.5 35.0]
        @test AxisKeys.dimnames(grouped2, 1) == :time
        @test AxisKeys.dimnames(grouped2, 2) == :component

        by_component = tmean(ka2, 1; dim = :component)
        @test by_component == [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0]

        @testset "output coordinate" begin
            @test axiskeys(grouped2, 1) == Millisecond.([0, 2])
            @test axiskeys(grouped2, 2) == [1, 2]
        end
    end
end
