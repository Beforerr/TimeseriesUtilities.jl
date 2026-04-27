@testitem "smooth" begin
    using Dates
    using DimensionalData
    using TimeseriesUtilities

    data = [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0]
    times = [0.0, 1.0, 3.0, 6.0]

    @test smooth(data, 3; dim = 1) == [1.5 15.0; 2.0 20.0; 3.0 30.0; 3.5 35.0]
    @test smooth(data, times, 2.0; dim = 1) == [1.0 10.0; 1.5 15.0; 3.0 30.0; 4.0 40.0]
    @test smooth(data', times, 2.0; dim = 2) == [1.0 1.5 3.0 4.0; 10.0 15.0 30.0 40.0]
    @test smooth(data', 3) == [1.5 2.0 3.0 3.5; 15.0 20.0 30.0 35.0]

    explicit_left = smooth(data, times, (2.0, 0.0); dim = 1)
    @test isnan(explicit_left[1, 1])
    @test isnan(explicit_left[1, 2])

    t0 = DateTime(2020)
    datetimes = t0 .+ Second.([0, 1, 4, 5])
    da = DimArray([1.0, 2.0, 4.0, 8.0], (Ti(datetimes),))
    smoothed = smooth(da, Second(2))

    @test smoothed == [1.0, 1.5, 4.0, 6.0]
    @test parent(dims(smoothed, Ti)) == datetimes
end
