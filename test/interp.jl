@testitem "tinterp interpolation" begin
    using Dates, DimensionalData
    using TimeseriesUtilities

    times = [DateTime(2020, 1, 1), DateTime(2020, 1, 2), DateTime(2020, 1, 3)]
    da = DimArray(0:2, (Ti(times),))
    @test tinterp(da, DateTime(2020, 1, 1, 12)) ≈ 0.5

    t2 = [DateTime(2020, 1, 1, 6), DateTime(2020, 1, 2, 18)]
    res2 = tinterp(da, t2)
    @test dims(res2, Ti).val == t2
    @test DimensionalData.span(dims(res2, Ti)) isa DimensionalData.Irregular
    t3 = DateTime(2020, 1, 1):Hour(6):DateTime(2020, 1, 2)
    @test DimensionalData.span(dims(tinterp(da, t3), Ti)) == DimensionalData.Regular(Hour(6))

    da3 = DimArray([1.0 4.0; 2.0 5.0; 3.0 6.0], (Ti(times), Y([10, 20])))
    @test tinterp(da3, t2) ≈ [1.25 4.25; 2.75 5.75]

    # a scalar time selects a sample
    t1 = DateTime(2020, 1, 1, 12)
    @test @inferred(tinterp(da3, t1)) ≈ [1.5, 4.5] && dims(tinterp(da3, t1)) == dims(da3, (Y,))
end

@testitem "DataInterpolations compatibility" begin
    using Dates, DimensionalData
    using DataInterpolations: LinearInterpolation, ExtrapolationType, CubicSpline
    using TimeseriesUtilities

    times = [DateTime(2020, 1, 1), DateTime(2020, 1, 2), DateTime(2020, 1, 3)]
    da = DimArray(0:2, (Ti(times),))

    dd = DimArray(0:2, (Ti(Date.(times)),))
    for interp in (LinearInterpolation, TimeseriesUtilities.LinearInterpolation)
        @test tinterp(dd, [DateTime(2020, 1, 2, 12)]; interp) == [1.5]
    end
    @test parent(tinterp(Timed([0.0, 1.0, 2.0], [0, 1, 2]), [0.5]; interp = LinearInterpolation)) == [0.5]

    before = DateTime(2019, 12, 31)
    @test tinterp(da, before; extrapolation = true) ≈ -1.0
    @test tinterp(da, before; interp = LinearInterpolation, extrapolation = ExtrapolationType.Linear) ≈ -1.0

    # per-channel fallback; a spline through linear data is exact
    t3 = [DateTime(2020, 1, i) for i in 1:5]
    da3 = DimArray([Float64(i * j) for i in 1:5, j in 1:3], (Ti(t3), Y(1:3)))
    @test parent(tinterp(da3, [DateTime(2020, 1, 2, 12)]; interp = CubicSpline)) ≈ [2.5 5.0 7.5]
end

@testitem "linear search hint" begin
    using TimeseriesUtilities: tinterp, LinearInterpolation
    # Repeated, backward and endpoint queries after a sorted run; nonlinear data so a
    # wrong segment changes the value. Reference: scalar path, which never uses the hint.
    t = cumsum(rand(50) .+ 0.1)
    u = sin.(t)
    q = [t[1]; sort(t[1] .+ (t[end] - t[1]) .* rand(200)); t[end]; t[3]; t[3]; t[20]]
    @test parent(tinterp(Timed(u, t), q)) ≈ LinearInterpolation(u, t).(q)
end

@testitem "FastInterpolations extension" begin
    using Dates, DimensionalData
    using FastInterpolations
    using TimeseriesUtilities

    t = collect(0.0:20.0)
    q = [0.5, 3.3, 19.9]
    A = randn(2, length(t), 3)  # time on dim 2, channels on dims 1 and 3
    for f in (cubic_interp, pchip_interp)  # `Series` path and per-channel path
        ref = [f(t, A[p, :, r], x) for p in 1:2, x in q, r in 1:3]
        @test parent(tinterp(Timed(A, t; dim = 2), q; interp = f)) ≈ ref
        @test parent(tinterp(Timed(permutedims(A, (2, 1, 3)), t; dim = 1), q; interp = f)) ≈ permutedims(ref, (2, 1, 3))
    end

    @test parent(tinterp(Timed(t, t), [21.0]; interp = linear_interp, extrapolation = true)) ≈ [21.0]

    times = DateTime(2020) .+ Hour.(0:3)
    da = DimArray([0.0, 1.0, 4.0, 9.0], (Ti(times),))
    @test tinterp(da, DateTime(2020, 1, 1, 1, 30); interp = cubic_interp) ≈ cubic_interp(0:3, parent(da), 1.5)
end

@testitem "AxisKeys tinterp" begin
    using AxisKeys
    using Dates
    using TimeseriesUtilities

    times = [DateTime(2020, 1, 1), DateTime(2020, 1, 2), DateTime(2020, 1, 3)]
    t2 = [DateTime(2020, 1, 1, 6), DateTime(2020, 1, 2, 18)]
    ka2 = KeyedArray(
        [1.0 4.0; 2.0 5.0; 3.0 6.0];
        time = times,
        component = [10, 20],
    )
    res2 = tinterp(ka2, t2)
    @test res2 ≈ [1.25 4.25; 2.75 5.75]
    @test axiskeys(res2) == (t2, [10, 20])
    @test AxisKeys.dimnames(res2) == (:time, :component)

    resampled = tresample(KeyedArray(0.0:2.0; time = times), Hour(12))
    @test axiskeys(resampled, 1) == DateTime(2020, 1, 1):Hour(12):DateTime(2020, 1, 3)
    @test resampled ≈ 0.0:0.5:2.0
end

@testitem "tsync" begin
    using Dates, DimensionalData
    include("./setup.jl")

    da1, da2, da3 = workload_interp_setup()
    a_sync, b_sync, c_sync = tsync(da1, da2, da3)

    # common span is the intersection of all inputs; da3 is interpolated onto da1's times
    tc = [DateTime(2020, 1, 2), DateTime(2020, 1, 3)]
    @test parent(dims(a_sync, Ti)) == parent(dims(b_sync, Ti)) == parent(dims(c_sync, Ti)) == tc
    @test a_sync == [2, 3]
    @test b_sync == [10, 11]
    @test parent(c_sync) ≈ [5.5 9; 6.5 11]

    using JET
    @test_opt tsync(da1, da2, da3)
    @test_call tsync(da1, da2, da3)
end
