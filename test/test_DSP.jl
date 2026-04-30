@testitem "tfilter" begin
    using DimensionalData
    using Unitful
    using DSP

    fs = 100u"Hz"
    n = 512
    t = Ti((0:(n - 1)) ./ ustrip(fs) .* u"s")
    # 5 Hz signal well inside a 1–20 Hz bandpass
    signal = sin.(2π * 5 .* (0:(n - 1)) ./ ustrip(fs))
    da = DimArray(signal, (t,))

    result = tfilter(da, 1u"Hz", 20u"Hz")
    @test size(result) == size(da)
    @test dims(result, Ti) == dims(da, Ti)
    # passband signal should be mostly preserved (within 5% after transient edges)
    interior = (n ÷ 10):(9n ÷ 10)
    @test maximum(abs.(result[interior] .- signal[interior])) < 0.05
end