using BenchmarkTools
import DataInterpolations
using DimensionalData
using TimeseriesUtilities
using TimeseriesUtilities: tinterp, tsync
include("../test/setup.jl")

const SUITE = BenchmarkGroup()

const da_bench = DimArray(rand(1000, 3), (Ti(1:1000), Y(1:3)))
const t_bench = rand(1:1000, 32)

SUITE["tinterp"] = BenchmarkGroup()
SUITE["tinterp"]["internal_linear"] = @benchmarkable tinterp($da_bench, $t_bench; interp = TimeseriesUtilities.LinearInterpolation)
SUITE["tinterp"]["data_interpolations_linear"] = @benchmarkable tinterp($da_bench, $t_bench; interp = DataInterpolations.LinearInterpolation)

const interp_nans_dimarray_bench = let data = rand(1000, 3)
    data[250:250:750, :] .= NaN
    DimArray(data, (Ti(1:1000), Y(1:3)))
end

SUITE["tinterp_nans"] = @benchmarkable tinterp_nans($interp_nans_dimarray_bench)

const da1_bench, da2_bench, da3_bench = workload_interp_setup(128)

SUITE["tsync"] = @benchmarkable tsync($da1_bench, $da2_bench, $da3_bench)

const smooth_data_bench = rand(1000, 3)
const smooth_times_bench = cumsum(rand(1000))
const smooth_dimarray_bench = DimArray(smooth_data_bench, (Ti(smooth_times_bench), Y(1:3)))

SUITE["smooth"] = BenchmarkGroup()
SUITE["smooth"]["sample_window"] = @benchmarkable smooth($smooth_data_bench, 25; dim = 1)
SUITE["smooth"]["time_window"] = @benchmarkable smooth($smooth_data_bench, $smooth_times_bench, 0.1; dim = 1)
SUITE["smooth"]["dimarray_time_window"] = @benchmarkable smooth($smooth_dimarray_bench, 0.1)
