@testitem "IntervalRange and tsplit" begin
    using Dates
    using Chairmarks

    # Test tsplit with integer n
    result = tsplit(1.0, 10, 3)
    @test collect(result) == [(1.0, 4.0), (4.0, 7.0), (7.0, 10.0)]

    result = tsplit(Date("2020-02-29"), Date("2022-03-01"), 4)
    @test result == [
        (Date("2020-02-29"), Date("2020-08-30")),
        (Date("2020-08-30"), Date("2021-03-01")),
        (Date("2021-03-01"), Date("2021-08-31")),
        (Date("2021-08-31"), Date("2022-03-01")),
    ]

    # Test tsplit with step size dt
    result = tsplit(0.0, 10.0, 2.0)
    @test length(result) == 5
    @test result[1] == (0.0, 2.0)
    @test result[5] == (8.0, 10.0)

    # Test tsplit with DateTime and Month period
    result = tsplit(DateTime("2021-01-01"), DateTime("2021-03-01"), Month)
    @test result == [
        (DateTime("2021-01-01"), DateTime("2021-02-01")),
        (DateTime("2021-02-01"), DateTime("2021-03-01")),
    ]

    eom = IntervalRange(Date("2021-01-31"), Date("2021-04-01"), Month(1))
    @test collect(eom) == [
        (Date("2021-01-31"), Date("2021-02-28")),
        (Date("2021-02-28"), Date("2021-03-31")),
        (Date("2021-03-31"), Date("2021-04-01")),
    ]
    @test eom[2] == (Date("2021-02-28"), Date("2021-03-31"))
    @test (@b collect(eom)).allocs <= 2

    quarter_split = IntervalRange(Date("2021-01-31"), Date("2021-08-01"), Quarter(1))
    @test collect(quarter_split) == [
        (Date("2021-01-31"), Date("2021-04-30")),
        (Date("2021-04-30"), Date("2021-07-31")),
        (Date("2021-07-31"), Date("2021-08-01")),
    ]

    year_split = IntervalRange(Date("2020-02-29"), Date("2022-03-01"), Year(1))
    @test year_split[3] == (Date("2022-02-28"), Date("2022-03-01"))
    @test_throws BoundsError year_split[4]

    # Test tsplit with DateTime and Day period
    result = tsplit(DateTime("2021-01-01"), DateTime("2021-01-05"), Day)
    @test result == [
        (DateTime("2021-01-01"), DateTime("2021-01-02")),
        (DateTime("2021-01-02"), DateTime("2021-01-03")),
        (DateTime("2021-01-03"), DateTime("2021-01-04")),
        (DateTime("2021-01-04"), DateTime("2021-01-05")),
    ]

    # Test tsplit with DateTime and Hour period
    result = tsplit(DateTime("2021-01-01T00:00:00"), DateTime("2021-01-01T03:00:00"), Hour)
    @test result == [
        (DateTime("2021-01-01T00:00:00"), DateTime("2021-01-01T01:00:00")),
        (DateTime("2021-01-01T01:00:00"), DateTime("2021-01-01T02:00:00")),
        (DateTime("2021-01-01T02:00:00"), DateTime("2021-01-01T03:00:00")),
    ]

    using JET
    @test_call tsplit(1.0, 10.0, 3)
    @test_call tsplit(0.0, 10.0, 2.0)
    @test_call tsplit(DateTime("2021-01-01"), DateTime("2021-03-01"), Month)
end
