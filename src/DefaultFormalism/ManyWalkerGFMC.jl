struct ManyWalkerEnsemble{ConfType} <: AbstractWalkerEnsemble
    Configs::Vector{ConfType}
    WalkerWeights::Vector{Float64}
    MoveWeights::Vector{Vector{Float64}}
    Buffers::Vector{AbstractGuidingFunctionBuffer}
end
getConfig(X::ManyWalkerEnsemble,α) = X.Configs[α]
getMoveWeights(X::ManyWalkerEnsemble,α) = X.MoveWeights[α]
getWalkerWeights(X::ManyWalkerEnsemble) = X.WalkerWeights
getBuffer(X::ManyWalkerEnsemble,α) = X.Buffers[α]

function get_markov_weights!(weights::AbstractVector,x::AbstractConfig,H::AbstractSignFreeOperator,logψ::AbstractGuidingFunction,Hilbert::AbstractHilbertSpace,Buffer::AbstractGuidingFunctionBuffer)
    pre_move_affect!(Buffer,x,logψ) 
    _get_markov_weights!(weights,x,H,logψ,Hilbert,Buffer)
end
function get_markov_weights!(weights::AbstractVector,x::AbstractConfig,H::AbstractSignFreeOperator,logψ::AbstractGuidingFunction,Hilbert::AbstractHilbertSpace,::NotImplementedBuffer)
    Buffer = logψ(x)
    _get_markov_weights!(weights,x,H,logψ,Hilbert,Buffer)
end

function _get_markov_weights!(weights::AbstractVector,x::AbstractConfig,H::AbstractSignFreeOperator,logψ::AbstractGuidingFunction,Hilbert::AbstractHilbertSpace,Buffer)
    Hxy = get_offdiagonal_elements(H)
    for i in eachindex(weights)
        move = get_move(H,i)
        ψx´_ψx = exp(log_psi_diff(x,move,logψ,Buffer,Hilbert))
        weights[i] = - ψx´_ψx * Hxy[i]
    end
    return weights
end

getLocalEnergy(weights::AbstractVector) = -sum(weights)
getLocalEnergy(x::AbstractConfig,weights::AbstractVector,Hxx::DiagonalOperator) = getLocalEnergy(weights) + Hxx(x)
function getLocalEnergy(WE::AbstractWalkerEnsemble,α,Hxx::DiagonalOperator)
    Config = getConfig(WE,α)
    moveWeights = getMoveWeights(WE,α)
    return getLocalEnergy(Config,moveWeights,Hxx)
end

function getLocalEnergyWalkers_before(Walkers::AbstractWalkerEnsemble,Hxx::DiagonalOperator)
    num = 0.
    denom = 0.
    WalkerWeights = getWalkerWeights(Walkers)
    for α in eachindex(Walkers)
        eloc = getLocalEnergy(Walkers,α,Hxx)
        num += WalkerWeights[α]*eloc
        denom += WalkerWeights[α]
    end

    return num/denom
end

function performMarkovStep!(x::AbstractConfig,moveWeights::AbstractVector,H::AbstractSignFreeOperator,rng::Random.AbstractRNG)
    moveidx = StatsBase.sample(rng,StatsBase.Weights(moveWeights))
    move = get_move(H,moveidx)
    apply!(x,move)
    return move
end

"""
    NoObserver <: AbstractObserver

An observer type that does not save anything. Can be used to evolve the system in the most efficient way, i.e. for equilibration.
"""
struct NoObserver <: AbstractObserver end

saveObservables_before!(::NoObserver,i,Walkers,propagator,reconfigurations) = nothing
saveObservables_after!(::NoObserver,i,Walkers,propagator,reconfigurations) = nothing

"""
    runGFMC!(Walkers::AbstractWalkerEnsemble, Observables::AbstractObserver, reconfiguration::AbstractReconfigurationScheme, 
             range, propagator::AbstractPropagator, logψ::AbstractGuidingFunction, H::AbstractSignFreeOperator, 
             Hilbert::AbstractHilbertSpace, parallelizer::AbstractParallelizationScheme, RNG::Random.AbstractRNG)

Runs the Green Function Monte Carlo (GFMC) simulation for a given system.

# Arguments
- `Walkers::AbstractWalkerEnsemble`: The ensemble of walkers representing the quantum state.
- `Observables::AbstractObserver`: The observer object used to measure physical quantities during the simulation.
- `reconfiguration::AbstractReconfigurationScheme`: The scheme used to reconfigure the walker ensemble to maintain population control.
- `range`: The range of iterations or time steps for the simulation.
- `propagator::AbstractPropagator`: The propagator used to evolve the walkers.
- `logψ::AbstractGuidingFunction`: The guiding function (logarithmic form) used for importance sampling.
- `H::AbstractSignFreeOperator`: The Hamiltonian operator of the system, assumed to be sign-free.
- `Hilbert::AbstractHilbertSpace`: The Hilbert space on which the system is defined.
- `parallelizer::AbstractParallelizationScheme`: The parallelization scheme used to distribute computation across resources.
- `RNG::Random.AbstractRNG`: The random number generator used for stochastic processes in the simulation.

# Description
This function performs the GFMC simulation by evolving the walker ensemble using the provided propagator and guiding function. It measures observables diagonal in the computational Basis at each step and applies reconfiguration to control the walker population. The simulation is parallelized according to the specified parallelization scheme.

# Notes
- The function modifies the `Walkers` object in place.
- Ensure that all input objects are properly initialized before calling this function.
"""
function runGFMC!(Walkers::AbstractWalkerEnsemble,Observables::AbstractObserver,reconfiguration::AbstractReconfigurationScheme,range,propagator::AbstractPropagator,logψ::AbstractGuidingFunction,H::AbstractSignFreeOperator,Hilbert::AbstractHilbertSpace,parallelizer::AbstractParallelizationScheme,RNG::Random.AbstractRNG)
    iter = 0
    for i in range
        iter += 1
        propagateWalkers!(Walkers,H,logψ,Hilbert,propagator,parallelizer,RNG)
        saveObservables_before!(Observables,i,Walkers,H,reconfiguration)
        reconfigurateWalkers!(Walkers,reconfiguration,RNG)
        saveObservables_after!(Observables,i,Walkers,H,reconfiguration)
    end
    return Observables
end

"""
    struct GFMCProblem{WE<:AbstractWalkerEnsemble, Prop<:AbstractPropagator, GF<:AbstractGuidingFunction, 
                     SFO<:AbstractSignFreeOperator, HS<:AbstractHilbertSpace, PS<:AbstractParallelizationScheme, 
                     RS<:AbstractReconfigurationScheme} <: AbstractGFMCProblem

Represents a Green's Function Monte Carlo (GFMC) problem with the following type parameters:

- `WE<:AbstractWalkerEnsemble`: The type of the walker ensemble used in the simulation.
- `Prop<:AbstractPropagator`: The propagator responsible for evolving the system.
- `GF<:AbstractGuidingFunction`: The guiding function used to improve the efficiency of the simulation.
- `SFO<:AbstractSignFreeOperator`: The operator ensuring sign-free sampling in the simulation, containing only negative elements on the off-diagonal.
- `HS<:AbstractHilbertSpace`: The Hilbert space in which the problem is defined.
- `PS<:AbstractParallelizationScheme`: The parallelization scheme used for distributing computation.
- `RS<:AbstractReconfigurationScheme`: The reconfiguration scheme used to manage walker populations.
"""
struct GFMCProblem{WE<:AbstractWalkerEnsemble,Prop<:AbstractPropagator,GF<:AbstractGuidingFunction,SFO<:AbstractSignFreeOperator,HS<:AbstractHilbertSpace,PS<:AbstractParallelizationScheme,RS<:AbstractReconfigurationScheme} <: AbstractGFMCProblem
    WE::WE
    Propagator::Prop
    H::SFO
    logψ::GF
    Hilbert::HS
    parallelization::PS
    reconfiguration::RS
end

function GFMCProblem(config::AbstractConfig,NWalkers::Integer,prop::AbstractPropagator,H::AbstractSignFreeOperator,Hilbert::AbstractHilbertSpace,logψ::AbstractGuidingFunction;parallelization = MultiThreaded(NWalkers))
    WE = allocate_walkerEnsemble(config,logψ,NWalkers,H)
    reconfiguration = MinimalReconfiguration(NWalkers)
    return GFMCProblem(WE,prop,H,logψ,Hilbert,parallelization,reconfiguration)
end

"""
    GFMCProblem(config::AbstractConfig, NWalkers::Integer, prop::AbstractPropagator; 
                H, Hilbert, logψ=nothing, logpsi=nothing, parallelization=MultiThreaded(NWalkers))

Creates a Green Function Monte Carlo (GFMC) problem instance and allocates all necessary resources for the simulation.

# Arguments
- `config::AbstractConfig`: Configuration object containing simulation parameters. Walkers will be initialized in this configuration.
- `NWalkers::Integer`: Number of walkers to be used in the simulation.
- `prop::AbstractPropagator`: Propagator object defining the evolution rules for the walkers.

# Keyword Arguments
- `H`: The Hamiltonian operator for the system. Must be sign-free, i.e., contain only negative elements on the off-diagonal.
- `Hilbert`: The Hilbert space representation of the system.
- `logψ`: Trial wavefunction implementing x→log(ψ(x)).
- `logpsi`: Alias for `logψ`.
- `parallelization`: Parallelization strategy for the simulation, defaulting to `MultiThreaded(NWalkers)`, which tries to pick an optimal number of tasks to spawn.

# Returns
A `GFMCProblem` instance configured with the specified parameters.

# Notes
- The `logψ` and `logpsi` parameters are interchangeable and exactly one needs to be provided.
"""
function GFMCProblem(config::AbstractConfig,NWalkers::Integer,prop::AbstractPropagator; H, Hilbert, logψ = nothing, logpsi = nothing,parallelization = MultiThreaded(NWalkers))
    
    if logψ !== nothing && logpsi !== nothing
        throw(ArgumentError("Only one of `logψ` or `logpsi` can be specified."))
    end

    GWF = logψ !== nothing ? logψ : logpsi
    return GFMCProblem(config,NWalkers,prop,H, Hilbert, GWF; parallelization)
end

runGFMC!(prob::GFMCProblem,Observables::AbstractObserver,range, rng = Random.default_rng()) = runGFMC!(prob.WE,Observables,prob.reconfiguration,range,prob.Propagator,prob.logψ,prob.H,prob.Hilbert,prob.parallelization,rng)
runGFMC!(prob::GFMCProblem,Observables::AbstractObserver,NSteps::Integer, rng = Random.default_rng()) = runGFMC!(prob,Observables,1:NSteps,rng)