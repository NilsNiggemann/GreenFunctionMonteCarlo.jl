"""
    ContinuousTimePropagator{T<:AbstractFloat} <: AbstractPropagator

A structure representing a imaginary time projection of a trial wavefunction, also referred to as the "Continuous time limit". 

# Type Parameters
- `T<:AbstractFloat`: The floating-point type used for numerical computations 
  (e.g., `Float64`, `Float32`).

# Supertype
- `AbstractPropagator`: This structure is a subtype of `AbstractPropagator`, 
  indicating that it implements the required interface for propagators in the 
  Green Function Monte Carlo framework.
# Fields 
- `dτ::T`: The time step for the propagation, which is a floating-point value. A small value may be inefficient in exploring the Hilbert space, while a large value will lead to a more unstable propagation. A good starting point is `dτ = 0.1`.
- `w_avg_estimate::T`: An estimate of the average weight to reduce floating point errors. Ideally given by the exact ground state energy of the system.
"""
struct ContinuousTimePropagator{T<:AbstractFloat} <: AbstractPropagator 
    dτ::T
    w_avg_estimate::T
end

struct InfinitePropagationTimeError <: Exception
    walker_index::Int
    τ_step::Float64
    el_x::Float64
    H_xx::Float64
    τleft::Float64
    max_move_weight::Float64
end

# InfinitePropagationTimeError(walker_index, τ_step, el_x, H_xx, τleft, max_move_weight) = InfinitePropagationTimeError(walker_index, Float64(τ_step), Float64(el_x), Float64(H_xx), Float64(τleft), Float64(max_move_weight))

function Base.showerror(io::IO, err::InfinitePropagationTimeError)
    print(io,
        "Infinite propagation time encountered for walker ", err.walker_index,
        ". Check guiding wavefunction/buffers. ",
        "(τ_step=", err.τ_step,
        ", el_x=", err.el_x,
        ", H_xx=", err.H_xx,
        ", τleft=", err.τleft,
        ", max_move_weight=", err.max_move_weight,
        ")"
    )
end

ContinuousTimePropagator(dτ::Real,w_avg_estimate::Real) = ContinuousTimePropagator(float(dτ),float(w_avg_estimate))
ContinuousTimePropagator(dτ::Real;w_avg_estimate=0.) = ContinuousTimePropagator(dτ,w_avg_estimate)

@inline propagateWalkers!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, propagator::ContinuousTimePropagator, parallelization::AbstractParallelizationScheme, RNGs::Vector{<:Random.AbstractRNG}) = continuos_time_propagation!(WE, H, logψ, Hilbert, propagator.dτ,propagator.w_avg_estimate, parallelization, RNGs)

"""
    continuos_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, dτ::Real, parallelization::MultiThreaded, RNG::Random.AbstractRNG = Random.default_rng())

Perform continuous time propagation on a walker ensemble for a fixed time step `dτ`.

# Arguments
- `WE::AbstractWalkerEnsemble`: The ensemble of walkers to be propagated.
- `H::AbstractSignFreeOperator`: The Hamiltonian operator used for propagation.
- `logψ::AbstractGuidingFunction`: The guiding function for the propagation.
- `Hilbert::AbstractHilbertSpace`: The Hilbert space in which the propagation occurs.
- `dτ::Real`: The time step for the propagation.
- `w_avg_estimate::Real`: An estimate of the average weight.
- `parallelization::MultiThreaded`: Parallelization settings for the propagation.
- `RNG::Random.AbstractRNG`: The random number generator to be used (default is `Random.default_rng()`).
"""
function continuos_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, dτ::Real, w_avg_estimate::Real, parallelization::MultiThreaded, RNGs::Vector{<:Random.AbstractRNG})
    
    batches = ChunkSplitters.chunks(eachindex(WE), n = num_tasks(parallelization),split= ChunkSplitters.RoundRobin())

    @sync for (i_chunk, αinds) in enumerate(batches)
        rng = RNGs[i_chunk]
         
        Threads.@spawn for α in αinds
            continuos_time_propagation_walker!(WE, α, H, logψ, Hilbert, dτ, w_avg_estimate, rng)
        end
    end
end

function continuos_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, dτ::Real, w_avg_estimate::Real, parallelization::BatchMultiThreaded, RNGs::Vector{<:Random.AbstractRNG})
    αinds = eachindex(WE)
    Nw = length(αinds)
    minbatch = clamp(Nw ÷ num_tasks(parallelization), 1, Nw)
    Polyester.@batch minbatch=minbatch per=thread for α in αinds
        task_idx = polyester_get_task_local_idx(α, minbatch)
        continuos_time_propagation_walker!(WE, α, H, logψ, Hilbert, dτ, w_avg_estimate, RNGs[task_idx])
    end
end

function continuos_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, dτ::Real, w_avg_estimate::Real, parallelization::SingleThreaded, RNGs::Vector{<:Random.AbstractRNG})
    for α in eachindex(WE)
        continuos_time_propagation_walker!(WE, α, H, logψ, Hilbert, dτ, w_avg_estimate, first(RNGs))
    end
end

function continuos_time_propagation_walker!(WE::AbstractWalkerEnsemble, α::Int, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, dτ::Real, w_avg_estimate::Real,RNG::Random.AbstractRNG)
    Config = getConfig(WE, α)
    GWFBuffer = getBuffer(WE, α)
    moveWeights = getMoveWeights(WE, α)
    log_w = 0.0
    get_markov_weights!(moveWeights, Config, H, logψ, Hilbert, GWFBuffer)
    Hxx = get_diagonal(H)
    H_xx = Hxx(Config)
    el_x = H_xx + getLocalEnergy(moveWeights)
    τleft = dτ
    while τleft > 0
        ξ = rand(RNG)
        τ_step = min(τleft, log(1 - ξ) / (el_x - H_xx))
        τleft -= τ_step
        if isinf(τleft)
            throw(InfinitePropagationTimeError(α,τ_step,el_x,H_xx,τleft,maximum(moveWeights)))
        end
        log_w += -τ_step * el_x
        if τleft > 0
            last_move = performMarkovStep!(Config, moveWeights, H, RNG)
            post_move_affect!(GWFBuffer, Config, last_move, logψ)
            get_markov_weights!(moveWeights, Config, H, logψ, Hilbert, GWFBuffer)

            H_xx = Hxx(Config)
            el_x = H_xx + getLocalEnergy(moveWeights)
        end
    end
    w = exp(log_w - dτ * w_avg_estimate)
    WalkerWeights = getWalkerWeights(WE)
    WalkerWeights[α] = w
    localEnergies = getLocalEnergies(WE)
    localEnergies[α] = el_x
    return nothing
end