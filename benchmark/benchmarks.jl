# Benchmark
using Chairmarks
using DimensionalData
using TimeseriesUtilities: tinterp, tsync, workload_interp_setup

da_bench = DimArray(rand(1000, 3), (Ti(1:1000), Y(1:3)))
t_bench = rand(1:1000, 32)
@info "tinterp" @b(tinterp(da_bench, t_bench))
# 6.458 μs (42 allocs: 2.250 KiB)

da1, da2, da3 = workload_interp_setup(128)
@info "tsync" @b(tsync(da1, da2, da3))
# 10.542 μs (78 allocs: 7.844 KiB)
