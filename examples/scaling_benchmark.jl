"""
Parallel-scaling benchmark for GreenFunctionMonteCarlo.jl
==========================================================

This script measures wall-clock time for a fixed GFMC simulation as the
number of tasks (and therefore threads) increases, and compares it against
the ideal linear speedup.

Run with multiple threads to see the benefit, e.g.:
    julia --threads 8 examples/scaling_benchmark.jl

Key conditions for near-perfect scaling
-----------------------------------------
1. **No allocations in the hot loop** – the `sample_fast` helper avoids the
   heap allocation that `StatsBase.Weights(...)` incurs on every Markov step.
2. **Small working set per walker** – choose N_sites small enough that the
   per-walker buffers (Jastrow h_i vector, v_ij matrix) fit in per-core L2
   cache so threads do not compete for memory bandwidth.
3. **Many walkers per thread** – use NWalkers ≫ nthreads so the load
   imbalance (random while-loop iterations per walker) averages out.
4. **Default TaskLocalRNG** – Julia's `TaskLocalRNG` is zero-overhead for
   threading: every spawned task gets its own copy of the RNG state, so there
   is no lock contention on the RNG object itself.
"""

using GreenFunctionMonteCarlo
using GreenFunctionMonteCarlo.LinearAlgebra
using Random

# ---------------------------------------------------------------------------
# Model: transverse-field Ising chain  H = -h Σ σ^x_i + J Σ σ^z_i σ^z_{i+1}
# At the critical point h = J = 1.
# Using HardCoreConstraint (spin-1/2 ≡ bosons with occupancy 0 or 1).
# ---------------------------------------------------------------------------

const N_SITES   = 16          # small enough that h_i fits comfortably in L1
const N_WALKERS = Threads.nthreads() * 200   # ≫ nthreads → good load balance
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

    # trivial (flat) guiding function → no Jastrow buffer cost
    logψ = EqualWeightSuperposition()

    config = BosonConfig(Hilbert)

    # energy estimate from a quick VMC run (reduces weight variance)
    prop = ContinuousTimePropagator(0.1; w_avg_estimate = -1.0)

    prob = GFMCProblem(config, NWalkers, prop; logψ, H, Hilbert, parallelization)
    return prob
end

# ---------------------------------------------------------------------------
# Benchmark helper: measure median wall time over a few repetitions
# ---------------------------------------------------------------------------

function bench_run(prob; nreps = 3)
    times = Float64[]
    for _ in 1:nreps
        t = @elapsed runGFMC!(prob, NoObserver(), N_BENCH; logger = nothing)
        push!(times, t)
    end
    return minimum(times)   # use minimum to reduce scheduling noise
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
    prob_st = build_problem(N_SITES, N_WALKERS; parallelization = SingleThreaded())
    runGFMC!(prob_st, NoObserver(), N_WARMUP; logger = nothing)  # equilibrate
    t_single = bench_run(prob_st)
    @printf("  SingleThreaded:  %.3f s\n\n", t_single)

    # Multi-threaded with increasing task counts
    println("Benchmarking MultiThreaded…")
    println("  tasks  time(s)  speedup  ideal")
    println("  -----  -------  -------  -----")

    for ntasks in unique([1, 2, 4, 8, nthreads, 2nthreads, 4nthreads])
        ntasks > 4 * N_WALKERS && continue   # nonsensical
        prob_mt = build_problem(N_SITES, N_WALKERS;
                                parallelization = MultiThreaded(ntasks))
        runGFMC!(prob_mt, NoObserver(), N_WARMUP; logger = nothing)
        t = bench_run(prob_mt)
        speedup = t_single / t
        ideal   = min(Float64(ntasks), Float64(nthreads))
        @printf("  %5d  %7.3f  %7.2fx  %4.1fx\n", ntasks, t, speedup, ideal)
    end

    println()
    println("""
Tips for better scaling
-----------------------
• Run with more threads: `julia --threads auto examples/scaling_benchmark.jl`
• Increase N_WALKERS relative to nthreads for better load balancing.
• Keep N_SITES small so per-walker buffers stay in CPU cache.
• Use the default TaskLocalRNG (default_rng()) — it has zero lock overhead
  because each spawned task gets its own independent RNG state.
• Avoid passing a shared non-task-local RNG (e.g. MersenneTwister) to the
  parallel path: those objects require a lock and serialize all rand() calls.
""")
end

# Printf is in stdlib; import it
import Printf: @printf

main()
