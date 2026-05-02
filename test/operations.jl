@testitem "tsort, tclip, tview" begin
    using DimensionalData

    times = [3.0, 1.0, 2.0]
    t = Ti(times)
    y = Y(1:3)
    da = rand(t, y)
    result = tsort(da)
    @test parent(dims(result, Ti)) == sort(times)
    @test result[Ti(1)] == da[Ti(2)]
    @test result[Ti(3)] == da[Ti(1)]

    # tclip
    @test_throws ArgumentError tclip(da, 1.0, 2.0) == da[2:3, :]
    @test tclip(result, 1.0, 2.0) == parent(da[2:3, :])

    # tview
    @test_throws ArgumentError tview(da, 1.0, 2.0)
    @test tview(result, 1.0, 2.0) == parent(da[2:3, :])

    ds = DimStack(
        (
            a = DimArray([30.0, 10.0, 20.0], (t,)),
            b = DimArray([30.0 300.0 3000.0; 10.0 100.0 1000.0; 20.0 200.0 2000.0], (t, y)),
        )
    )
    sorted_ds = tsort(ds)
    @test dims(sorted_ds, Ti).val == sort(times)
    @test sorted_ds.a == [10.0, 20.0, 30.0]
    @test_throws ArgumentError tclip(ds, 1.0, 2.0)
    clipped_ds = tclip(sorted_ds, 1.0, 2.0)
    @test dims(clipped_ds, Ti).val == [1.0, 2.0]
    @test clipped_ds.a == [10.0, 20.0]
    @test clipped_ds.b == [10.0 100.0 1000.0; 20.0 200.0 2000.0]
    @test tview(sorted_ds, 1.0, 2.0).a == clipped_ds.a

    # Benchmark
    using Chairmarks
    verbose = false
    tclip_bench = @b tclip($result, 1.0, 2.0)
    tview_bench = @b tview($result, 1.0, 2.0)
    @test tclip_bench.allocs > tview_bench.allocs
    @test tclip_bench.time > tview_bench.time
    verbose && @info tclip_bench, tview_bench

    using JET
    @test_opt tsort(da)
    @test_call tsort(da)
    @test_opt tclip(result, 1.0, 2.0)
    @test_call tclip(result, 1.0, 2.0)
    @test_opt tview(result, 1.0, 2.0)
    @test_call tview(result, 1.0, 2.0)
end

@testitem "dropna" begin
    using DimensionalData

    t = Ti(2:5)
    y = Y([:x, :y])
    ds = DimStack(
        (
            a = DimArray([1.0, NaN, 3.0, 4.0], (t,)),
            b = DimArray([1.0 5.0; 2.0 NaN; 3.0 7.0; 4.0 8.0], (t, y)),
        )
    )

    a_clean = dropna(ds.a)
    @test a_clean == [1.0, 3.0, 4.0]
    @test a_clean.dims[1].val == [2, 4, 5]
    result = dropna(ds)
    @test dims(result, Ti).val == [2, 4, 5]
    @test result.a == [1.0, 3.0, 4.0]
    @test result.b == [1.0 5.0; 3.0 7.0; 4.0 8.0]
end

@testitem "AxisKeys time operations" begin
    using AxisKeys
    using Chairmarks

    times = [3.0, 1.0, 2.0, 5.0, 4.0]
    values = [30.0, 10.0, 20.0, 50.0, 40.0]
    ka = KeyedArray(values; time = times)

    sorted = tsort(ka)
    @test axiskeys(sorted, 1) == sort(times)
    @test sorted == [10.0, 20.0, 30.0, 40.0, 50.0]

    @test_throws ArgumentError tclip(ka, 2.0, 4.0)
    clipped = tclip(sorted, 2.0, 4.0)
    @test clipped == [20.0, 30.0, 40.0]
    @test axiskeys(clipped, 1) == [2.0, 3.0, 4.0]

    viewed = tview(sorted, 2.0, 4.0)
    @test viewed == clipped
    @test (@b tview($sorted, 2.0, 4.0)).allocs ≤ 1
    @test axiskeys(viewed, 1) == axiskeys(clipped, 1)

    @test tselect(sorted, 3.4) == 30.0
    @test tselect(ka, 3.4) == 30.0
    @test tselect(sorted, 3.4, 0.5) == 30.0
    @test ismissing(tselect(sorted, 1.5, 0.25))

    masked = tmask(sorted, 2.0, 4.0)
    @test isequal(masked, [10.0, NaN, NaN, NaN, 50.0])
    @test axiskeys(masked, 1) == axiskeys(sorted, 1)

    ka2 = KeyedArray(reshape(1.0:15.0, 5, 3); time = [1.0, 2.0, 3.0, 4.0, 5.0], space = [:x, :y, :z])
    @test tview(ka2, 2.0, 4.0) == ka2[2:4, :]
    @test (@b tview($ka2, 2.0, 4.0)).allocs ≤ 1
    @test tselect(ka2, 2.2) == ka2[2, :]
    @test axiskeys(tclip(ka2, 2.0, 4.0), 1) == [2.0, 3.0, 4.0]
    @test axiskeys(tclip(ka2, 2.0, 4.0), 2) == [:x, :y, :z]

    ka3 = KeyedArray([100.0, 200.0, 300.0]; time = [3.0, 4.0, 5.0])
    c1, c2 = tclips(sorted, ka3)
    @test axiskeys(c1, 1) == axiskeys(c2, 1) == [3.0, 4.0, 5.0]
    v1, v2 = tviews(sorted, ka3)
    @test v1 == c1
    @test v2 == c2
end


@testitem "tclips, tviews" begin
    using DimensionalData

    da1 = DimArray([10.0, 20.0, 30.0, 40.0, 50.0], (Ti([1.0, 2.0, 3.0, 4.0, 5.0]),))
    da2 = DimArray([1.0, 2.0, 3.0, 4.0, 5.0], (Ti([2.0, 3.0, 4.0, 5.0, 6.0]),))
    da3 = DimArray([100.0], (Ti([11.0]),))

    # common range is [2.0, 5.0]
    c1, c2 = tclips(da1, da2)
    @test dims(c1, Ti).val == [2.0, 3.0, 4.0, 5.0]
    @test dims(c2, Ti).val == [2.0, 3.0, 4.0, 5.0]
    @test c1 == [20.0, 30.0, 40.0, 50.0]

    # explicit trange overrides
    c1e, c2e = tclips(da1, da2; trange = (3.0, 4.0))
    @test dims(c1e, Ti).val == [3.0, 4.0]

    # tviews — same logic, view semantics
    v1, v2 = tviews(da1, da2)
    @test dims(v1, Ti).val == [2.0, 3.0, 4.0, 5.0]
    @test dims(v2, Ti).val == [2.0, 3.0, 4.0, 5.0]

    @test_throws ArgumentError tclips(da1, da3)

    using JET
    @test_opt tclips(da1, da2)
    @test_call tclips(da1, da2)
end


@testitem "tselect" begin
    using DimensionalData

    # Test with numeric time dimension
    times = [1.0, 3.0, 5.0, 7.0, 9.0]

    @test_throws MethodError tsort(times)
    @test_throws MethodError tclip(times, 3.0, 7.0)
    @test_throws MethodError tview(times, 3.0, 7.0)
    @test_throws MethodError tselect(times, 5.0)

    values = [10.0, 20.0, 30.0, 40.0, 50.0]
    da = DimArray(values, (Ti(times),))

    # Test exact match
    @test tselect(da, 5.9) == tselect(da, 5.0, 0.5) == tselect(da, 5.2, 0.5) == tselect(da, 4.8, 0.5) == 30.0
    # Test closest match at edge of range
    @test tselect(da, 5.5, 0.5) == tselect(da, 4.5, 0.5) == 30.0

    # Test no match within range
    @test ismissing(tselect(da, 2.0, 0.5))
    @test ismissing(tselect(da, 8.0, 0.5))

    using JET
    @test_opt tselect(da, 5.0, 0.5)
    @test_call tselect(da, 5.0, 0.5)
end

@testitem "tmask!, tmask" begin
    using DimensionalData

    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    da = DimArray([10.0, 20.0, 30.0, 40.0, 50.0], (Ti(times),))

    # tmask! (mutating)
    da_mut = copy(da)
    tmask!(da_mut, 2.0, 4.0)
    @test isequal(da_mut, [10.0, NaN, NaN, NaN, 50.0])
    # original unchanged
    @test da[Ti(At(2.0))] == 20.0

    # tmask (non-mutating copy)
    masked = tmask(da, 2.0, 4.0)
    @test da[Ti(At(2.0))] == 20.0
    @test isequal(masked, [10.0, NaN, NaN, NaN, 50.0])

    # array of intervals
    da_multi = tmask(da, [(1.0, 2.0), (4.0, 5.0)])
    @test isequal(da_multi, [NaN, NaN, 30.0, NaN, NaN])

    using JET
    @test_opt tmask!(copy(da), 2.0, 4.0)
    @test_call tmask!(copy(da), 2.0, 4.0)
end

@testitem "tshift" begin
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    values = [10.0, 20.0, 30.0, 40.0, 50.0]

    @testset "DimensionalData" begin
        using DimensionalData
        da = DimArray(values, (Ti(times),))

        # Test shift by 2.0
        result = tshift(da, 2.0)
        @test parent(dims(result, Ti)) == [-1.0, 0.0, 1.0, 2.0, 3.0]

        # Test default shift with 2D array
        y = Y(1:3)
        da2d = rand(Ti(times), y)
        result2d = tshift(da2d)
        @test parent(dims(result2d, Ti)) == times .- 1.0
    end

    @testset "AxisKeys" begin
        using AxisKeys
        ka = KeyedArray(values; time = times)

        # Test shift by 2.0
        result = tshift(ka, 2.0)
        @test axiskeys(result, 1) == [-1.0, 0.0, 1.0, 2.0, 3.0]

        # Test with 2D array
        ka2d = KeyedArray(rand(5, 3); time = times, space = 1:3)
        result2d = tshift(ka2d)
        @test axiskeys(result2d, 1) == times .- 1.0
        @test axiskeys(result2d, 2) == 1:3
    end
end
