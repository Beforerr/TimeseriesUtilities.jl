# Benchmark extrema/median implementations for resolution computation.
# Run manually: julia --project=test test/bench_resolution.jl
using Dates
using Chairmarks
using Test

t = Millisecond.(0:3:100000)

# cadence
using TimeseriesUtilities: resolution
using SpaceDataModel: cadence
@info "Cadence vs Resolution" (@b cadence($t), resolution($t))
@assert cadence(t) == resolution(t)

# median and extrema benchmarks
dt = diff(t)

import VectorizedStatistics
import TimeseriesUtilities.NaNStatistics, Statistics

for f in (:extrema, :median)
    @info "Benchmarking $f implementations:"
    f2_sym = Symbol("v$f")
    f3_sym = Symbol("nan$f")
    @eval begin
        f1 = ts -> Statistics.$f(reinterpret(Int, ts))
        f2 = ts -> VectorizedStatistics.$f2_sym(reinterpret(Int, ts))
        f3 = ts -> NaNStatistics.$f3_sym(reinterpret(Int, ts))
    end
    b1, b2, b3 = @b f1(dt), f2(dt), f3(dt)
    @test f1(dt) == f2(dt) == f3(dt)
    times = (b1.time, b2.time, b3.time)
    names = ("Statistics.$f", "VectorizedStatistics.$f2_sym", "NaNStatistics.$f3_sym")
    fastest_idx = argmin(times)
    @info "Times: $(names .=> times)"
    @info "Fastest $f implementation: $(names[fastest_idx])"
end
