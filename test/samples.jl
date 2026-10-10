@testitem "mapsamples" begin
    using Dates
    using DimensionalData
    using AxisKeys: KeyedArray, axiskeys, dimnames

    A = [1.0 4.0; 2.0 5.0; 3.0 6.0]
    ts = [DateTime(2020), DateTime(2021)]
    shift(x, t) = x .+ year(t)

    @test mapsamples(shift, A, ts; dim = 2) == A .+ [2020 2021]
    @test mapsamples(shift, permutedims(A), ts; dim = 1) == permutedims(A .+ [2020 2021])
    @test mapsamples(shift, A, ts[1]; dim = 2) == A .+ 2020
    @test mapsamples(shift, A, Tuple(ts); dim = 2) == A .+ [2020 2021]
    @test_throws DimensionMismatch mapsamples(shift, A, ts[1:1]; dim = 2)
    @test_throws ArgumentError mapsamples(shift, A, ts[1]; dim = 3)
    @test collect(eachsample(A, ts; dim = 2)) == [(A[:, 1], ts[1]), (A[:, 2], ts[2])]
    @test_throws DimensionMismatch mapsamples((x, t) -> x[1:t], A, [2, 1]; dim = 2)
    # element types widen like `map`
    r = mapsamples((x, t) -> t == 1 ? missing : sum(x), A, 1:2; dim = 2)
    @test r isa Vector{Union{Missing, Float64}} && isequal(r, [missing, 15.0])
    r = mapsamples((x, t) -> t == 1 ? Int.(x) : x / 2, A, 1:2; dim = 2)
    @test eltype(r) == Real && r == [1 2; 2 2.5; 3 3]
    @test mapsamples(shift, zeros(3, 0), DateTime[]; dim = 2) == zeros(3, 0)
    @test mapsamples((x, t) -> sum(x), zeros(3, 0), DateTime[]; dim = 2) isa Vector{Float64}
    # result with fewer axes than a sample: time goes last
    B = rand(3, 4, 5)
    @test mapsamples((x, t) -> x .+ t, B, 1:4; dim = 2) == B .+ reshape(1:4, 1, 4)
    @test mapsamples((x, t) -> vec(sum(x; dims = 2)), B, 1:5; dim = 3) == dropdims(sum(B; dims = 2); dims = 2)
    @test mapsamples((x, t) -> vec(sum(x; dims = 2)), zeros(3, 4, 0), 1:0; dim = 3) == zeros(0, 0)

    da = DimArray(permutedims(A), (Ti(ts), Y([:x, :y, :z])))
    r = mapsamples(shift, da)
    @test DimensionalData.dims(r) == DimensionalData.dims(da) && parent(r) == permutedims(A .+ [2020 2021])
    r = mapsamples((x, t) -> sum(x), da)
    @test r isa DimVector && lookup(r, Ti) == ts
    @test lookup(mapsamples((x, t) -> sum(x), set(da, Ti => Dim{:time})), Dim{:time}) == ts
    r = mapsamples((x, t) -> x[1:2], da)
    @test lookup(r, Ti) == ts && size(r) == (2, 2) && DimensionalData.dims(r, 2) isa DimensionalData.Dimensions.AnonDim
    r = mapsamples((x, t) -> sum(x), da[1:0, :])
    @test r isa DimVector && isempty(r) && eltype(lookup(r, Ti)) == DateTime

    ka = KeyedArray(A; c = [:x, :y, :z], time = ts)
    r = mapsamples(shift, ka)
    @test dimnames(r) == (:c, :time) && axiskeys(r) == axiskeys(ka)
    r = mapsamples((x, t) -> year(t), ka)
    @test dimnames(r) == (:time,) && r == [2020, 2021]
    r = mapsamples((x, t) -> x[1:2], ka)
    @test dimnames(r) == (:_, :time) && axiskeys(r, :time) == ts
end

@testitem "SpaceDataModel variables" begin
    using Dates
    import SpaceDataModel as SDM

    struct TimedVar{T, N, A <: AbstractArray{T, N}, D} <: SDM.AbstractDataVariable{T, N}
        data::A
        dims::D
    end
    SDM.tdimnum(::TimedVar) = 1

    ts = DateTime(2020) .+ Day.(0:2)
    v = TimedVar([1.0 2.0; 3.0 4.0; 5.0 6.0], (ts, [:a, :b]))
    @test times(v) == ts
    @test tclip(v, ts[2], ts[3]) == parent(v)[2:3, :]
    @test mapsamples((x, t) -> sum(x) + day(t), v) == [3.0 + 1, 7.0 + 2, 11.0 + 3]
end

@testitem "mapsamples performance" begin
    using Dates, Chairmarks, DimensionalData

    scalar(x, t) = x[1] * x[2] + Dates.value(t)
    vector(x, t) = x
    times_(n) = DateTime(2020) .+ Second.(1:n)
    data(n, d) = d == 2 ? rand(3, n) : rand(n, 3)
    # a runtime `dim` union-splits for any rank
    @inferred mapsamples(vector, rand(2, 3, 4), 1:3; dim = Ref(2)[])
    for f in (scalar, vector), d in (1, 2)
        @inferred mapsamples(f, data(4, d), times_(4); dim = d)
        allocs(n) = (A = data(n, d); ts = times_(n); (@b mapsamples($f, $A, $ts; dim = $d)).allocs)
        @test allocs(1000) == allocs(2000)
    end
    shared(n) = (A = data(n, 2); (@b mapsamples(scalar, $A, DateTime(2020); dim = 2)).allocs)
    @test shared(1000) == shared(2000)
    # wrapped arrays: the time dim comes from the dimension types and must stay static
    da = DimArray(rand(3, 4, 5), (X(1:3), Y(1:4), Ti(times_(5))))
    @inferred mapsamples((x, t) -> sum(x), da)
end
