using BenchmarkTools
import DataInterpolations
using DimensionalData
using Random
using StaticArrays: SVector
using TimeseriesUtilities
using TimeseriesUtilities: tinterp, tsync
include("../test/setup.jl")

Random.seed!(42)

const SUITE = BenchmarkGroup()

const da_bench = DimArray(rand(1000, 3), (Ti(1:1000), Y(1:3)))
const t_bench = rand(1:1000, 32)

SUITE["tinterp"] = BenchmarkGroup()
SUITE["tinterp"]["internal_linear"] = @benchmarkable tinterp($da_bench, $t_bench; interp = TimeseriesUtilities.LinearInterpolation)
SUITE["tinterp"]["internal_linear_sorted_dense"] = @benchmarkable tinterp($da_bench, $(1.5:0.25:999.5))
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
SUITE["smooth"]["time_window"] = @benchmarkable smooth($(Timed(smooth_data_bench, smooth_times_bench; dim = 1)), 0.1)
SUITE["smooth"]["dimarray_time_window"] = @benchmarkable smooth($smooth_dimarray_bench, 0.1)

const samples_times_bench = collect(1.0:10_000.0)
sample_scalar(x, t) = sum(x) * t
sample_vector(x, t) = x .* t
sample_static(x, t) = SVector{3}(x) .* t

SUITE["mapsamples"] = BenchmarkGroup()
SUITE["mapsamples"]["scalar_time_last"] = @benchmarkable mapsamples(sample_scalar, $(rand(3, 10_000)), $samples_times_bench; dim = 2)
SUITE["mapsamples"]["vector_time_first"] = @benchmarkable mapsamples(sample_vector, $(rand(10_000, 3)), $samples_times_bench; dim = 1)
SUITE["mapsamples"]["static_time_first"] = @benchmarkable mapsamples(sample_static, $(rand(10_000, 3)), $samples_times_bench; dim = 1)
SUITE["mapsamples"]["static_time_last"] = @benchmarkable mapsamples(sample_static, $(rand(3, 10_000)), $samples_times_bench; dim = 2)
SUITE["mapsamples"]["dimarray_scalar"] = @benchmarkable mapsamples(sample_scalar, $(DimArray(rand(10_000, 3), (Ti(samples_times_bench), Y(1:3)))), $samples_times_bench)
