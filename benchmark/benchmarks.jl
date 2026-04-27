using BenchmarkTools
using DimensionalData
using TimeseriesUtilities
using TimeseriesUtilities: tinterp, tsync, workload_interp_setup

const SUITE = BenchmarkGroup()

const da_bench = DimArray(rand(1000, 3), (Ti(1:1000), Y(1:3)))
const t_bench = rand(1:1000, 32)

SUITE["tinterp"] = @benchmarkable tinterp($da_bench, $t_bench)

const da1_bench, da2_bench, da3_bench = workload_interp_setup(128)

SUITE["tsync"] = @benchmarkable tsync($da1_bench, $da2_bench, $da3_bench)

const smooth_data_bench = rand(1000, 3)
const smooth_times_bench = cumsum(rand(1000))
const smooth_dimarray_bench = DimArray(smooth_data_bench, (Ti(smooth_times_bench), Y(1:3)))

SUITE["smooth"] = BenchmarkGroup()
SUITE["smooth"]["sample_window"] = @benchmarkable smooth($smooth_data_bench, 25; dim = 1)
SUITE["smooth"]["time_window"] = @benchmarkable smooth($smooth_data_bench, $smooth_times_bench, 0.1; dim = 1)
SUITE["smooth"]["dimarray_time_window"] = @benchmarkable smooth($smooth_dimarray_bench, 0.1)
