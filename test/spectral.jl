@testitem "tfilter raw arrays" begin
    using DSP

    fs = 100
    n = 512
    t = (0:(n - 1)) ./ fs
    signal = sin.(2π * 5 .* t)
    interior = (n ÷ 10):(9n ÷ 10)

    result = parent(tfilter(Timed(signal, t), 1, 20))
    @test maximum(abs.(result[interior] .- signal[interior])) < 0.05

    multi = permutedims(hcat(signal, cos.(2π * 5 .* t)))
    multi_result = parent(tfilter(Timed(multi, t; dim = 2), 1, 20))
    @test maximum(abs.(multi_result[1, interior] .- signal[interior])) < 0.05
end

@testitem "tfilter axis containers" begin
    using AxisKeys
    using DimensionalData
    using DSP

    fs = 100
    n = 512
    t = (0:(n - 1)) ./ fs
    signal = sin.(2π * 5 .* t)
    interior = (n ÷ 10):(9n ÷ 10)

    da = DimArray(signal, (Ti(t),))
    result = tfilter(da, 1, 20)

    keyed = KeyedArray(signal; time = t)
    keyed_result = tfilter(keyed, 1, 20)

    @test size(result) == size(keyed_result) == size(signal)
    @test dims(result, Ti) == dims(da, Ti)
    @test axiskeys(keyed_result, 1) == axiskeys(keyed, 1)
    @test maximum(abs.(result[interior] .- signal[interior])) < 0.05

    multi = DimArray(permutedims(hcat(signal, cos.(2π * 5 .* t))), (Dim{:comp}([:x, :y]), Ti(t)))
    multi_result = tfilter(multi, 1, 20; dim = Ti)
    @test size(multi_result) == size(multi)
    @test dims(multi_result, Ti) == dims(multi, Ti)
    @test dims(multi_result, Dim{:comp}) == dims(multi, Dim{:comp})
end

@testitem "pspectrum raw arrays" begin
    using Dates
    using DSP

    ts = DateTime(2020):Second(1):(DateTime(2020) + Second(99))
    signal = sin.(1:100)
    raw = pspectrum(Timed(signal, ts); nfft = 16)

    @test size(parent(raw)) == (9, 11)
    @test times(raw)[1] == DateTime("2020-01-01T00:00:08")
    @test TimeseriesUtilities.axiskeys(raw, 1)[1] == 0.0

    ts_time = Time(0):Nanosecond(500_000_000):(Time(0) + Nanosecond(500_000_000 * 99))
    @test times(pspectrum(Timed(signal, ts_time); nfft = 16))[1] == Time(0) + Second(4)

    multi = permutedims(hcat(signal, cos.(1:100)))
    multi_spec = pspectrum(Timed(multi, ts; dim = 2); nfft = 16)
    @test size(parent(multi_spec)) == (9, 11, 2)
    @test multi_spec.keys[1:2] == raw.keys
end

@testitem "pspectrum axis containers" begin
    using AxisKeys
    using Dates
    using DimensionalData
    using DSP

    ts = DateTime(2020):Second(1):(DateTime(2020) + Second(99))
    da = DimArray(sin.(1:100), (Ti(ts),))

    y = pspectrum(da; nfft = 16)
    raw = pspectrum(Timed(parent(da), ts); nfft = 16)
    y_default = pspectrum(da)
    keyed = KeyedArray(sin.(1:100); time = ts)
    keyed_spec = pspectrum(keyed; nfft = 16)

    @test size(y) == size(keyed_spec) == (9, 11)
    @test parent(y) == parent(keyed_spec) == parent(raw)
    @test dims(y, Dim{:frequency}).val == axiskeys(keyed_spec, :frequency) == raw.keys[1]
    @test dims(y, Ti).val == axiskeys(keyed_spec, :time) == raw.keys[2]
    @test size(y_default, 1) == 129
    @test dims(y, Ti)[1] == DateTime("2020-01-01T00:00:08")
    @test dims(y, Dim{:frequency})[1] == 0.0

    multi = DimArray(hcat(sin.(1:100), cos.(1:100)), (Ti(ts), Dim{:comp}([:x, :y]; metadata = Dict(:units => "nT"))))
    multi_spec = pspectrum(multi; nfft = 16)
    @test DimensionalData.metadata(dims(multi_spec, Dim{:comp}))[:units] == "nT"
    multi_t2 = DimArray(permutedims(parent(multi)), (Dim{:comp}([:x, :y]), Ti(ts)))
    multi_t2_spec = pspectrum(multi_t2; dim = Ti, nfft = 16)
    keyed_multi = KeyedArray(permutedims(parent(multi)); comp = [:x, :y], time = ts)
    keyed_multi_spec = pspectrum(keyed_multi; dim = :time, nfft = 16)

    @test size(multi_spec) == size(multi_t2_spec) == size(keyed_multi_spec) == (9, 11, 2)
    @test parent(multi_spec) == parent(keyed_multi_spec)
    @test dims(multi_spec, Dim{:frequency}).val == axiskeys(keyed_multi_spec, :frequency) == raw.keys[1]
    @test dims(multi_spec, Ti).val == axiskeys(keyed_multi_spec, :time) == raw.keys[2]
    @test dims(multi_spec, Dim{:comp}).val == dims(multi_t2_spec, Dim{:comp}).val == axiskeys(keyed_multi_spec, :comp) == [:x, :y]
    @test AxisKeys.dimnames(keyed_multi_spec) == (:frequency, :time, :comp)
end
