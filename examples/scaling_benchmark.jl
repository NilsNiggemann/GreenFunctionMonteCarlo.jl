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

const N_SITES   = 16          # small enough that buffers fit comfortably in L1
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
    Hxx = DiagOperator(x -> sum(σz(x[i]) * σz(x[i+1]) for i in 1:N-1))
    H   = localOperator(moves, offdiag, Hxx, Hilbert)

    # trivial (flat) guiding function -> no Jastrow buffer cost
    logψ = EqualWeightSuperposition()

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
    print("Warming up (compiling + equilibrating)... ")
    prob_warmup = build_problem(N_SITES, N_WALKERS)
    runGFMC!(prob_warmup, NoObserver(), N_WARMUP; logger = nothing)
    println("done.\n")

    # Baseline: single-threaded
    print("Benchmarking SingleThreaded... ")
    prob_st = build_problem(N_SITES, N_WALKERS; parallelization = SingleThreaded())
    runGFMC!(prob_st, NoObserver(), N_WARMUP; logger = nothing)  # equilibrate
    t_single = bench_run(prob_st)
    @printf("%.3f s\n\n", t_single)

    # Multi-threaded with increasing task counts
    println("Benchmarking MultiThreaded...")
    println("  tasks  time(s)  speedup  ideal  Amdahl_max")
    println("  -----  -------  -------  -----  ----------")

    t_1task = nothing
    for ntasks in unique([1, 2, 4, 8, nthreads, 2nthreads, 4nthreads])
        ntasks > 4 * N_WALKERS && continue   # nonsensical
        prob_mt = build_problem(N_SITES, N_WALKERS;
                                parallelization = MultiThreaded(ntasks))
        runGFMC!(prob_mt, NoObserver(), N_WARMUP; logger = nothing)
        t = bench_run(prob_mt)
        speedup = t_single / t
        ideal   = min(Float64(ntasks), Float64(nthreads))
        if ntasks == 1
            t_1task = t
            @printf("  %5d  %7.3f  %7.2fx  %4.1fx        ---\n", ntasks, t, speedup, ideal)
        else
            # Amdahl fit from 1-task and nthreads-task measurements:
            # T(n) = t_par/n + t_seq
            # t_seq = (n*T(n) - T(1)) / (n - 1)   (only meaningful at large n)
            if ntasks <= nthreads && t_1task !== nothing
                n = Float64(ntasks)
                t_seq = (n * t - t_1task) / (n - 1)
                t_seq = max(t_seq, 0.0)  # clamp numerical noise
                S_max = t_1task / t_seq
                @printf("  %5d  %7.3f  %7.2fx  %4.1fx   %6.2fx\n", ntasks, t, speedup, ideal, S_max)
            else
                @printf("  %5d  %7.3f  %7.2fx  %4.1fx        ---\n", ntasks, t, speedup, ideal)
            end
        end
    end

    if t_1task !== nothing && nthreads >= 2
        # Estimate sequential fraction at n = nthreads
        n  = Float64(nthreads)
        # Use 1-task and nthreads-task timings for the best Amdahl estimate
        prob_nt = build_problem(N_SITES, N_WALKERS;
                                parallelization = MultiThreaded(nthreads))
        runGFMC!(prob_nt, NoObserver(), N_WARMUP; logger = nothing)
        t_nt = bench_run(prob_nt)
        t_seq = (n * t_nt - t_1task) / (n - 1)
        t_seq = max(t_seq, 0.0)
        f_seq = t_seq / t_1task   # sequential fraction
        S_max = 1 / f_seq

        println()
        @printf("Amdahl analysis (1 task vs %d tasks):\n", nthreads)
        @printf("  Estimated sequential fraction : %.1f%%\n", 100 * f_seq)
        @printf("  Theoretical max speedup       : %.2fx\n", S_max)
        println()
        println("The sequential fraction comes primarily from `reconfigurateWalkers!`.")
        println("This step does global resampling and cannot be parallelised.")
        println("The `minimizeReconfiguration!` helper has been optimised to use a")
        println("pre-allocated inverse buffer (no Dict allocation per step),")
        println("which reduces GC pressure that would otherwise stall all threads.")
    end

    println()
    println("""
Tips for better scaling
-----------------------
* Run with more threads: `julia --threads auto examples/scaling_benchmark.jl`
* Increase N_WALKERS relative to nthreads for better load balancing.
* Keep N_SITES small so per-walker buffers stay in CPU cache.
* Use the default TaskLocalRNG (default_rng()) -- it has zero lock overhead
  because each spawned task gets its own independent RNG state.
* Avoid passing a shared non-task-local RNG (e.g. MersenneTwister) to the
  parallel path: those objects require a lock and serialize all rand() calls.
* Consider increasing dτ to reduce the average number of inner-loop Markov
  steps per walker (fewer calls into performMarkovStep! per GFMC step).
""")
end

main()
