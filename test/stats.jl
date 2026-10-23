@testitem "timeseries statistics" begin
    using Dates
    using Statistics

    @testset "tmean" begin
        t0 = DateTime(2011)
        offsets = Millisecond.(0:9)
        ts = t0 .+ offsets
        @test tmean(ts) == t0 + Microsecond(4500)

        using Durations: Timestamp
        tn = Timestamp{Nanosecond}(t0) .+ Nanosecond.([0, 1, 2, 3, 10])
        @test tmean(tn) === tn[4]
        @test tmedian(tn) === tn[3]
        @test parent(tmedian(Timed([ts ts], 1:10; dim = 1))) == [ts[5], ts[5]]
        @test parent(tmean(Timed(ts, 1:10), 5)) == [ts[3], ts[7], ts[10]]  # windows [0, 5), [5, 10), [10, 15)
    end

    @testset "DimensionalData" begin
        using DimensionalData

        t = Ti(Millisecond.(0:3))
        y = Y(1:2)
        da1 = rand(t)
        da2 = rand(t, y)

        @test tmean(da1) ≈ mean(da1)
        @test tmean(da1, Millisecond(2)) == [mean(da1[1:2]), mean(da1[3:4])]
        @test @inferred(tmean(da2)) ≈ vec(mean(da2, dims = 1))
        @test tmean(da2, Millisecond(2)) == [mean(parent(da2)[1:2, :], dims = 1); mean(parent(da2)[3:4, :], dims = 1)]

        # tmedian
        @test tmedian(da1) == median(da1)
        @test tmedian(da1, Millisecond(2)) == [median(da1[1:2]), median(da1[3:4])]
        @test tmedian(da2) == vec(median(da2, dims = 1))

        @test tstd(da1) ≈ std(da1)

        # DimStack
        @test tmean(DimStack((da1, da2))) == (; layer1 = tmean(da1), layer2 = tmean(da2))

        @testset "output coordinate" begin
            grouped2 = tmean(da2, Millisecond(2))
            @test grouped2.dims[1].val == Millisecond.([0, 2])
            @test grouped2.dims[2].val == 1:2
        end
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

        @test @inferred(tmean(ka2)) == [2.5, 25.0] && axiskeys(tmean(ka2), 1) == [1, 2]

        by_component = tmean(ka2, 1; dim = :component)
        @test by_component == [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0]

        @testset "output coordinate" begin
            @test axiskeys(grouped2, 1) == Millisecond.([0, 2])
            @test axiskeys(grouped2, 2) == [1, 2]
        end
    end
end

@testitem "windowed tstat" begin
    using Dates

    x = Timed([1.0 10.0; 2.0 20.0; 3.0 30.0; 9.0 90.0], [0.0, 1.0, 2.5, 7.0]; dim = 1)
    # overlapping windows; empty ones stay so starts are regular
    r = tmean(x, 2.0; period = 3.0)
    @test times(r) == 0.0:2.0:6.0
    @test isequal(parent(r), [2.0 20.0; 3.0 30.0; NaN NaN; 9.0 90.0])
    @test times(tsum(x, 2.0; start_by = :datapoint)) == 0.0:2.0:6.0
    @test times(tsum(Timed([1.0, 2.0], [0.5, 3.5]), 2.0; start_by = :datapoint)) == 0.5:2.0:2.5
    @test_throws ArgumentError tsum(x, 2.0; start_by = :start)
    # `f` sees one channel and may return anything
    @test parent(tstat(extrema, x, 4.0))[:, 1] == [(1.0, 3.0), (9.0, 9.0)]
    @test tstat(length, x) == [4, 4]

    # channels land back in place when time is a middle dimension
    A = reshape(1.0:24, 2, 4, 3)
    w = tsum(Timed(A, [0.0, 1.0, 2.5, 7.0]; dim = 2), 2.0)
    @test parent(w) == cat(sum(A[:, 1:2, :]; dims = 2), A[:, 3:3, :], zeros(2, 1, 3), A[:, 4:4, :]; dims = 2)

    using DimensionalData
    ds = DimStack((a = DimArray([1.0, 2.0, 3.0], Ti([0.0, 1.0, 2.5])),))
    @test parent(tsum(ds, 2.0; period = 3.0).a) == [6.0, 3.0]

    t = DateTime(2020) .+ Minute.([1, 2, 5])
    @test times(tmedian(Timed([1.0, 2.0, 3.0], t), Minute(2))) == DateTime(2020):Minute(2):DateTime(2020, 1, 1, 0, 4)
end
