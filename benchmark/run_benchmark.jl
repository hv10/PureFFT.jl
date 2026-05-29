#!/usr/bin/env julia
# Benchmark script for PureFFT: compare fft (four radix/method combos) vs dft
# Produces CSV and PNG in scripts/output/

using Pkg
# Make sure the package source is discoverable for local dev
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))
using PureFFT
using BenchmarkTools
using CSV
using DataFrames
using StatsPlots
using Random
using FFTW
using FFTA

# Make sure the package source is discoverable for local dev
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

# reproducible RNG
Random.seed!(1234)

output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)

# A list of increasing highly-composite / highly-factored numbers.
# Adjust as desired; keep sizes small enough so naive DFT finishes in reasonable time.
# sizes = [60, 120, 240, 360, 720, 840, 1260, 1680, 2520]
sizes = [360]

# Four method/radix combos
combos = (
    (method=:dit, rad=:min),
    (method=:dif, rad=:min),
    (method=:dit, rad=:max),
    (method=:dif, rad=:max),
)

results = DataFrame(N=Int[], method=String[], radix=String[], impl=String[], time_s=Float64[])

# Number of benchmark samples per measurement (keeps runtime reasonable)
samples = 6

for N in sizes
    println("\n== N = $N ==")
    # create a random complex input
    x = randn(N) .+ im * randn(N)

    # measure DFT once per N (naive)
    println(" Warm-up / compile: dft")
    dft(x; NN=N, inverse=false, normalize=false)
    println(" Benchmarking dft (N=$N)")
    bd = @benchmark dft($x; NN=$N, inverse=false, normalize=false) samples = samples
    dft_time = median(bd).time / 1e9
    println("  dft: $(round(dft_time, sigdigits=4)) s (median)")
    push!(results, (N, "dft", "-", "PureFFT", dft_time))

    # our fft variants
    for (combo) in combos
        method = combo.method
        rad = combo.rad
        println(" Warm-up / compile: fft (method=$method, rad=$rad)")
        try
            PureFFT.fft(x; method=method, rad=rad)
        catch e
            @warn "fft($N) with method=$method rad=$rad failed during warmup: $e"
            continue
        end

        println(" Benchmarking fft (method=$method, rad=$rad, N=$N)")
        b = @benchmark PureFFT.fft($x; method=$method, rad=$rad) samples = samples
        fft_time = median(b).time / 1e9
        println("  fft: $(round(fft_time, sigdigits=4)) s (median)")
        push!(results, (N, string(method), string(rad), "PureFFT", fft_time))
    end

    # FFTW
    println(" Warm-up: FFTW.fft")
    try
        FFTW.fft(x)
        bfw = @benchmark FFTW.fft($x) samples = samples
        fftw_time = median(bfw).time / 1e9
        println("  FFTW: $(round(fftw_time, sigdigits=4)) s (median)")
        push!(results, (N, "fft", "-", "FFTW", fftw_time))
    catch e
        @warn "FFTW failed for N=$N: $e"
    end

    # FFTA
    try
        @eval using FFTA
        println(" Warm-up: FFTA.fft")
        FFTA.fft(x)
        bffta = @benchmark FFTA.fft($x) samples = samples
        ffta_time = median(bffta).time / 1e9
        println("  FFTA: $(round(ffta_time, sigdigits=4)) s (median)")
        push!(results, (N, "fft", "-", "FFTA", ffta_time))
    catch
        # ignore if FFTA not available
        @warn "FFTA failed for N=$N: $e"
    end
end

csv_out = joinpath(output_dir, "benchmark_results.csv")
CSV.write(csv_out, results)
println("Wrote CSV: $csv_out")

# Plot: runtime (s) vs size; plot dft baseline and the other impls
plt = plot(title="PureFFT: runtime vs size", xlabel="N", ylabel="time (s)", legend=:topleft)

# dft baseline (unique per N)
dft_per_N = filter(row -> row.impl == "PureFFT" && row.method == "dft", results)
sort!(dft_per_N, :N)
# plot!(plt, dft_per_N.N, dft_per_N.time_s, lw=3, ls=:dash, color=:black, label="dft (naive)")

# plot each implementation / method
for impl in unique(results.impl)
    for m in unique(results.method)
        sub = filter(row -> row.impl == impl && row.method == m && row.method != "dft", results)
        if nrow(sub) == 0
            continue
        end
        sort!(sub, :N)
        label = string(impl, " ", m)
        plot!(plt, sub.N, sub.time_s, marker=:o, label=label)
    end
end

png_out = joinpath(output_dir, "benchmark_plot.png")
savefig(plt, png_out)
println("Wrote plot: $png_out")

# Compute speedups vs dft per N and write CSV
speedups = DataFrame(N=Int[], impl=String[], method=String[], speedup=Float64[])
for row in eachrow(results)
    if row.method == "dft" && row.impl == "PureFFT"
        continue
    end
    d = dft_per_N[dft_per_N.N.==row.N, :time_s][1]
    push!(speedups, (row.N, row.impl, row.method, d / row.time_s))
end

speedup_csv = joinpath(output_dir, "benchmark_speedups.csv")
CSV.write(speedup_csv, speedups)
println("Wrote speedup CSV: $speedup_csv")

println("\nDone. Results in: $output_dir")
