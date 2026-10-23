@testitem "JET" tags = [:jet] begin
    using JET, Dates, DimensionalData
    include("setup.jl")

    data = Timed((1:5) .^ 2, 1:5)
    A = DimArray([i^2 for i in 1:5, j in 1:3], (Ti(1.0:5.0), Y(1:3)))
    @test_opt ignored_modules = (Base,) tderiv(data)
    @test_call ignored_modules = (Base,) tderiv(data)
    @test_opt tderiv(A)
    @test_call ignored_modules = (Base,) tderiv(A)

    @test_call tsplit(1.0, 10.0, 3)
    @test_call tsplit(0.0, 10.0, 2.0)
    @test_call tsplit(DateTime("2021-01-01"), DateTime("2021-03-01"), Month)

    shift(x, t) = x .+ year(t)
    ts = [DateTime(2020), DateTime(2021)]
    @test_call mapsamples(shift, [1.0 4.0; 2.0 5.0; 3.0 6.0], ts; dim = 2)
    @test_opt mapsamples(shift, [1.0 4.0; 2.0 5.0; 3.0 6.0], ts; dim = 2)

    da1, da2, da3 = workload_interp_setup()
    @test_opt tsync(da1, da2, da3)
    @test_call tsync(da1, da2, da3)

    unsorted = rand(Ti([3.0, 1.0, 2.0]), Y(1:3))
    sorted = tsort(unsorted)
    @test_opt tsort(unsorted)
    @test_call tsort(unsorted)
    @test_opt tclip(sorted, 1.0, 2.0)
    @test_call tclip(sorted, 1.0, 2.0)
    @test_opt tview(sorted, 1.0, 2.0)
    @test_call tview(sorted, 1.0, 2.0)

    da = DimArray([10.0, 20.0, 30.0, 40.0, 50.0], (Ti(1.0:5.0),))
    @test_opt tselect(da, 3.0, 0.5)
    @test_call tselect(da, 3.0, 0.5)
    @test_opt tmask!(copy(da), 2.0, 4.0)
    @test_call tmask!(copy(da), 2.0, 4.0)
end
