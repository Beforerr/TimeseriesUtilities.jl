
@testitem "tmin, tmax, timerange" begin
    using Chairmarks
    using Dates
    using DimensionalData

    t = Ti(1:10000)
    x = rand(t)
    @test tminimum(x) == minimum(t)
    @test tmaximum(x) == maximum(t)
    @test timerange(x) == extrema(t)

    verbose = false

    for T in (Int, Date, DateTime)
        ts = T.(collect(1:10000))
        @test timerange(ts) == extrema(ts)
        b1 = @b timerange($ts)
        b2 = @b extrema($ts)
        if b1.time < b2.time
            verbose && @info "Acceleration ratio: $(b2.time / b1.time)"
        else
            @info "Deceleration ratio: $(b1.time / b2.time)"
        end
    end
end