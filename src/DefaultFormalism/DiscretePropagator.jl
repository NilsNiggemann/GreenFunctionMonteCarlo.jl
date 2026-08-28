"""
    DiscretePropagator{T<:AbstractFloat} <: AbstractPropagator

A structure representing a fixed-generation-count, discrete-time imaginary time projection of a trial wavefunction.

Each call to `propagateWalkers!` advances every walker by `nBranch` elementary Markov sub-steps. At each sub-step, the
walker either stays at its current configuration `x` (a "self-loop", with weight `Λ - H_xx(x)`) or moves to a
configuration `y` reachable via an off-diagonal matrix element (with weight `ψ(y)/ψ(x) * (-H_xy)`), sampled proportionally
to these weights. This implements the standard discretization of the imaginary-time propagator `G = Λ·I - H`, which
requires `Λ` to be chosen large enough that `Λ - H_xx(x) ≥ 0` for every configuration `x` reachable during the
simulation, so that all sampling weights stay non-negative.

# Type Parameters
- `T<:AbstractFloat`: The floating-point type used for numerical computations (e.g., `Float64`, `Float32`).

# Supertype
- `AbstractPropagator`: This structure is a subtype of `AbstractPropagator`, indicating that it implements the
  required interface for propagators in the Green Function Monte Carlo framework.

# Fields
- `Λ::T`: The constant shift in `G = Λ·I - H`. Must satisfy `Λ ≥ H_xx(x)` for all configurations `x` visited during the
  simulation, so that the self-loop weight `Λ - H_xx(x)` stays non-negative.
- `nBranch::Int`: The number of Markov sub-steps performed per call to `propagateWalkers!`. Plays the analogous role to
  `dτ` for `ContinuousTimePropagator`: a larger `nBranch` gives a finer, more accurate discretization of the projection
  at increased computational cost.
- `w_avg_estimate::T`: An estimate of the average branching factor per sub-step, used to reduce floating point errors.
  Unlike `ContinuousTimePropagator` (where `w_avg_estimate` enters as an additive log-shift and defaults to `0.`), here
  it enters as a *divisor* (`bx = total / w_avg_estimate`), so it must be strictly positive; it defaults to `1.`.
"""
struct DiscretePropagator{T<:AbstractFloat} <: AbstractPropagator
    Λ::T
    nBranch::Int
    w_avg_estimate::T
end
DiscretePropagator(Λ::Real, nBranch::Integer, w_avg_estimate::Real) = DiscretePropagator(float(Λ), Int(nBranch), float(w_avg_estimate))
DiscretePropagator(Λ::Real, nBranch::Integer; w_avg_estimate=1.) = DiscretePropagator(Λ, nBranch, w_avg_estimate)

@inline propagateWalkers!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, propagator::DiscretePropagator, parallelization::AbstractParallelizationScheme, RNGs::Vector{<:Random.AbstractRNG}) =
    discrete_time_propagation!(WE, H, logψ, Hilbert, propagator.Λ, propagator.nBranch, propagator.w_avg_estimate, parallelization, RNGs)

"""
    discrete_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, Λ::Real, nBranch::Integer, w_avg_estimate::Real, parallelization::AbstractParallelizationScheme, RNGs::Vector{<:Random.AbstractRNG})

Perform `nBranch` discrete-time Markov sub-steps of branching random walk on a walker ensemble.

# Arguments
- `WE::AbstractWalkerEnsemble`: The ensemble of walkers to be propagated.
- `H::AbstractSignFreeOperator`: The Hamiltonian operator used for propagation.
- `logψ::AbstractGuidingFunction`: The guiding function for the propagation.
- `Hilbert::AbstractHilbertSpace`: The Hilbert space in which the propagation occurs.
- `Λ::Real`: The constant shift in `G = Λ·I - H`.
- `nBranch::Integer`: The number of Markov sub-steps to perform.
- `w_avg_estimate::Real`: An estimate of the average branching factor.
- `parallelization::AbstractParallelizationScheme`: Parallelization settings for the propagation.
- `RNGs::Vector{<:Random.AbstractRNG}`: A vector of random number generators, one per task (see `duplicate_rng`).
"""
function discrete_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, Λ::Real, nBranch::Integer, w_avg_estimate::Real, parallelization::MultiThreaded, RNGs::Vector{<:Random.AbstractRNG})
    batches = ChunkSplitters.chunks(eachindex(WE), n = num_tasks(parallelization), split= ChunkSplitters.RoundRobin())

    @sync for (i_chunk, αinds) in enumerate(batches)
        rng = RNGs[i_chunk]

        Threads.@spawn for α in αinds
            discrete_time_propagation_walker!(WE, α, H, logψ, Hilbert, Λ, nBranch, w_avg_estimate, rng)
        end
    end
end
function discrete_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, Λ::Real, nBranch::Integer, w_avg_estimate::Real, parallelization::BatchMultiThreaded, RNGs::Vector{<:Random.AbstractRNG})
    αinds = eachindex(WE)
    Nw = length(αinds)
    minbatch = clamp(Nw ÷ num_tasks(parallelization), 1, Nw)
    Polyester.@batch minbatch=minbatch per=thread for α in αinds
        task_idx = polyester_get_task_local_idx(α, minbatch)
        discrete_time_propagation_walker!(WE, α, H, logψ, Hilbert, Λ, nBranch, w_avg_estimate, RNGs[task_idx])
    end
end
function discrete_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, Λ::Real, nBranch::Integer, w_avg_estimate::Real, parallelization::SingleThreaded, RNGs::Vector{<:Random.AbstractRNG})
    for α in eachindex(WE)
        discrete_time_propagation_walker!(WE, α, H, logψ, Hilbert, Λ, nBranch, w_avg_estimate, first(RNGs))
    end
end

function discrete_time_propagation_walker!(WE::AbstractWalkerEnsemble, α::Int, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, Λ::Real, nBranch::Integer, w_avg_estimate::Real, RNG::Random.AbstractRNG)
    Config = getConfig(WE, α)
    GWFBuffer = getBuffer(WE, α)
    moveWeights = getMoveWeights(WE, α)
    Hxx = get_diagonal(H)

    w = 1.0
    for _ in 1:nBranch
        get_markov_weights!(moveWeights, Config, H, logψ, Hilbert, GWFBuffer)
        selfWeight = Λ - Hxx(Config)
        total = sum(moveWeights) + selfWeight
        w *= total / w_avg_estimate

        if rand(RNG) * total >= selfWeight
            last_move = performMarkovStep!(Config, moveWeights, H, RNG)
            post_move_affect!(GWFBuffer, Config, last_move, logψ)
        end
    end

    get_markov_weights!(moveWeights, Config, H, logψ, Hilbert, GWFBuffer)
    el_x = Hxx(Config) + getLocalEnergy(moveWeights)

    WalkerWeights = getWalkerWeights(WE)
    WalkerWeights[α] = w
    localEnergies = getLocalEnergies(WE)
    localEnergies[α] = el_x
    return nothing
end
