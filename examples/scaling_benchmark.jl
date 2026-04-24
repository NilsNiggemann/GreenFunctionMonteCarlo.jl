"""
Parallel-scaling benchmark for GreenFunctionMonteCarlo.jl
==========================================================

This script measures wall-clock time for a fixed GFMC simulation as the
number of tasks (and therefore threads) increases, and compares it against
the ideal linear speedup.

Run with multiple threads to see the benefit, e.g.:
    julia --threads 8 examples/scaling_benchmark.jl

Understanding the scaling limit (Amdahl's law)
-----------------------------------------------
Each GFMC step has two phases:

  1. Walker propagation  — fully parallel, time ∝ 1/nthreads
  2. Walker reconfiguration — inherently sequential (global resampling),
     time is independent of nthreads

Let t_par = time for the parallel phase (1 thread) and
    t_seq = time for the sequential phase.

Then: T(n) ≈ t_par/n + t_seq

The maximum achievable speedup is:
  S_max = (t_par + t_seq) / t_seq  = 1 + t_par/t_seq

If t_seq / (t_par + t_seq) = f (the sequential fraction), then:
  S_max = 1/f

The benchmark below measures t_par and t_seq from the 1-task and 8-task
timings and prints the implied maximum speedup.

Key conditions for near-ideal scaling of the parallel phase
------------------------------------------------------------
1. No allocations in the hot loop — `sample_fast` avoids the heap
   allocation that `StatsBase.Weights(...)` caused on every Markov step.
2. No allocations in reconfiguration — `minimizeReconfiguration!` now uses
   a pre-allocated inverse buffer instead of allocating a Dict{Int,Int} on
   every GFMC step.  The Dict allocation was the dominant source of GC
   pressure: 200 steps × ~50 KB/Dict = ~10 MB of short-lived objects that
   trigger stop-the-world GC pauses and stall all threads.
3. Small working set per walker — keep N_sites small so the per-walker
   buffers fit in per-core L2 cache.
4. Many walkers per thread — NWalkers ≫ nthreads so load imbalance
   (random while-loop iterations per walker) averages out.
5. Default TaskLocalRNG — Julia's `TaskLocalRNG` is zero-overhead:
   each spawned task gets its own independent RNG state via seeded
   copying, so there is no lock contention on the RNG object.
"""

using GreenFunctionMonteCarlo
using GreenFunctionMonteCarlo.LinearAlgebra
using Random
import Printf: @printf

# ---------------------------------------------------------------------------
# Model: transverse-field Ising chain  H = -h Σ σ^x_i + J Σ σ^z_i σ^z_{i+1}
# At the critical point h = J = 1.
# Using HardCoreConstraint (spin-1/2 = bosons with occupancy 0 or 1).
# ---------------------------------------------------------------------------

const N_SITES   = 300          # small enough that buffers fit comfortably in L1
const N_WALKERS = Threads.nthreads() * 200   # >> nthreads -> good load balance
const N_WARMUP  = 50          # steps to discard (equilibration)
const N_BENCH   = 200         # steps to time

σz(n::Bool) = 1 - 2n

function build_problem(N, NWalkers; parallelization = MultiThreaded(Threads.nthreads()))
    Hilbert = BosonHilbertSpace(N, HardCoreConstraint())

    # off-diagonal: single spin flip on site i
    moves = eachcol(Bool.(LinearAlgebra.I(N)))
    offdiag = -ones(N)                         # -h, h = 1

    # diagonal: Ising interaction J = 1
    function Hxx_func(config)
        sum(σz(config[i]) * σz(config[i+1]) for i in eachindex(config)[1:end-1])
    end
    function Hxx_func(config,E0,lastmove)
        j = lastmove.inds[1]
        # Only sites i-1 and i can change energy, so we can compute the difference from the old energy
        ΔE = 0

        s(i) = 1 <= i <= length(config) ? σz(config[i]) : 0
        ΔE -= 2 * s(j)* (s(j-1) + s(j+1))
        return ΔE + E0
    end
    Hxx = DiagOperator(Hxx_func)
    # Hxx = DiagOperator(x -> 0.)
    H   = localOperator(moves, offdiag, Hxx, Hilbert)

    # trivial (flat) guiding function -> no Jastrow buffer cost
    logψ = (Jastrow(zeros(Float32, N), zeros(Float32, N,N)))   # logψ = exp(Σ h_i n_i) with h_i = 0
    # logψ = EqualWeightSuperposition()

    config = BosonConfig(Hilbert)

    # energy estimate from a quick VMC run (reduces weight variance)
    prop = ContinuousTimePropagator(0.1; w_avg_estimate = -1.0)

    prob = GFMCProblem(config, NWalkers, prop; logψ, H, Hilbert, parallelization)
    return prob
end

# ---------------------------------------------------------------------------
# Benchmark helper: measure minimum wall time over a few repetitions
# ---------------------------------------------------------------------------

function bench_run(prob; nreps = 3)
    times = Float64[]
    for _ in 1:nreps
        t = @elapsed runGFMC!(prob, NoObserver(), N_BENCH; logger = nothing)
        push!(times, t)
    end
    return minimum(times)   # minimum reduces OS scheduling noise
end

# ---------------------------------------------------------------------------
# Main benchmark loop
# ---------------------------------------------------------------------------

function main()
    nthreads = Threads.nthreads()
    println("Julia threads available: $nthreads")
    println("N_sites    = $N_SITES")
    println("N_walkers  = $N_WALKERS  ($(N_WALKERS ÷ nthreads) per thread)")
    println("N_steps    = $N_BENCH  (after $N_WARMUP warmup steps)")
    println()

    # Warmup run (JIT compilation + equilibration)
    println("Warming up (compiling + equilibrating)…")
    prob_warmup = build_problem(N_SITES, N_WALKERS)
    runGFMC!(prob_warmup, NoObserver(), N_WARMUP; logger = nothing)
    println("Done.\n")

    # Baseline: single-threaded
    println("Benchmarking SingleThreaded…")
    prob_st = build_problem(N_SITES, N_WALKERS, SingleThreaded())
    runGFMC!(prob_st, NoObserver(), N_WARMUP; logger = nothing)  # equilibrate
    t_single = bench_run(prob_st)
    @printf("  SingleThreaded:  %.3f s\n\n", t_single)

    # Multi-threaded with increasing task counts
    println("Benchmarking MultiThreaded…")
    println("  tasks tasks/threads  time(s)  speedup  s vs M(1)  ideal  efficiency")
    println("  ----- -------------  -------  -------  ---------  -----  ----------")

    t_M1 = 0
    # for ntasks in unique([1, nthreads])
    for ntasks in sort!(unique([1,2,4,8,16, nthreads, 2nthreads,3nthreads,4nthreads, 8nthreads,32nthreads]))
        ntasks > 4 * N_WALKERS && continue   # nonsensical
        prob_mt = build_problem(N_SITES, N_WALKERS, BatchMultiThreaded(ntasks,N_WALKERS))
        runGFMC!(prob_mt, NoObserver(), N_WARMUP; logger = nothing)
        t = bench_run(prob_mt)
        if ntasks == 1
            t_M1 = t
        end
        speedup = t_single / t

        speedup_vs_M1 = t_M1 / t

        ideal   = min(Float64(ntasks), Float64(nthreads))
        efficiency = speedup / ideal * 100
        @printf("  %5d %13.3f  %7.3f  %6.2fx %9.2fx  %4.1fx  %8.1f%%\n", ntasks,ntasks/nthreads, t, speedup, speedup_vs_M1, ideal, efficiency)
    end

end

# Printf is in stdlib; import it
function ideal_scaling_2(N_res,N_iter,print=false)
    res = zeros(N_res)
    Threads.@threads for i in 1:N_res
    # for i in 1:N_res
        t = @elapsed begin
            resi = 0.0
            # for j in 1:N_iter
            GreenFunctionMonteCarlo.LoopVectorization.@turbo for j in 1:N_iter
                resi += sin(sqrt(Float64(j)))
            end
            res[i] = resi
        end
        print && println("Thread $(Threads.threadid()) is working on index $i: t = $(round(t, digits=3))s")
    end
    return res
end
pinnedcores = fast_core_cpuids()
if isinteractive()
    
    pinnedcores = 16:27
    # pinnedcores = 27:-1:16
    # pinnedcores = 12:-1:0
end
println("Pinning threads to CPU cores: $(join(pinnedcores, ","))")
# pinthreads(0:27; nthreads = Threads.nthreads(), threadpool = :all, force = true, warn = false)
pinthreads(pinnedcores; nthreads = Threads.nthreads(), threadpool = :all, force = true, warn = false)
##
main()
# benchmark_ideal_scaling()
# pinthreads(:cores)
##
# threadinfo(color=true)
# ideal_scaling_2(Threads.nthreads(), 1000)
# @time ideal_scaling_2(Threads.nthreads(), 300000000,true)
# benchmark_fast_core_scaling()
# isinteractive() || exit()
##
nthreads = Threads.nthreads()
ntasks = 1*nthreads
# prob_mt = build_problem(N_SITES, N_WALKERS, MultiThreaded(ntasks))
prob_mt = build_problem(N_SITES, N_WALKERS, BatchMultiThreaded(ntasks,N_WALKERS))
# prob_mt = build_problem(N_SITES, N_WALKERS, Batch)
runGFMC!(prob_mt, NoObserver(), N_WARMUP; logger = nothing)
# @profview t = bench_run(prob_mt)
@time runGFMC!(prob_mt, NoObserver(), 3N_BENCH; logger = nothing)
@time runGFMC!(prob_mt, NoObserver(), 3N_BENCH; logger = nothing)
# @time runGFMC!(prob_mt, NoObserver(), N_BENCH; logger = nothing)
# display(@benchmark runGFMC!(prob_mt, NoObserver(), 3N_BENCH; logger = nothing))
##
@profview runGFMC!(prob_mt, NoObserver(), 3N_BENCH; logger = nothing)
