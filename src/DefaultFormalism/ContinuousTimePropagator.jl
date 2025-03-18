struct ContinuousTimePropagator{T<:AbstractFloat} <: AbstractPropagator 
    dτ::T
    w_avg_estimate::T
end
ContinuousTimePropagator(dτ::Real,w_avg_estimate::Real) = ContinuousTimePropagator(float(dτ),float(w_avg_estimate))
ContinuousTimePropagator(dτ::Real;w_avg_estimate=0.) = ContinuousTimePropagator(dτ,w_avg_estimate)

@inline propagateWalkers!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, propagator::ContinuousTimePropagator, parallelization::AbstractParallelizationScheme, RNG::Random.AbstractRNG = Random.default_rng()) = continuos_time_propagation!(WE, H, logψ, Hilbert, propagator.dτ,propagator.w_avg_estimate, parallelization, RNG)

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
function continuos_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, dτ::Real, w_avg_estimate::Real, parallelization::MultiThreaded, RNG::Random.AbstractRNG = Random.default_rng())
    
    batches = ChunkSplitters.chunks(eachindex(WE), n = parallelization.nTasks,split= ChunkSplitters.RoundRobin())

    @sync for (i_chunk, αinds) in enumerate(batches)
        Threads.@spawn for α in αinds
            continuos_time_propagation_walker!(WE, α, H, logψ, Hilbert, dτ, w_avg_estimate, RNG)
        end
    end
end
function continuos_time_propagation!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, dτ::Real, w_avg_estimate::Real, parallelization::SingleThreaded, RNG::Random.AbstractRNG = Random.default_rng())
    for α in eachindex(WE)
        continuos_time_propagation_walker!(WE, α, H, logψ, Hilbert, dτ, w_avg_estimate, RNG)
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
            @info "" τ_step el_x H_xx τleft maximum(moveWeights)
            error("Infinite propagation time encountered. Check for too large values in guiding wavefunction or its Buffers!")
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
end
