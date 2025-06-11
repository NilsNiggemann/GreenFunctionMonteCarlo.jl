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
        weights[i] = log_psi_diff(x, move, logψ, Buffer, Hilbert)
    end
    LoopVectorization.@turbo for i in eachindex(Hxy,weights)
        weights[i] = -Hxy[i]*exp(weights[i])
    end
    return weights
end

getLocalEnergy(weights::AbstractVector) = -sum(weights)
getLocalEnergy(x::AbstractConfig,weights::AbstractVector,Hxx::DiagonalOperator) = getLocalEnergy(weights) + Hxx(x)
function getLocalEnergy(x::AbstractConfig,H::AbstractSignFreeOperator,logψ::AbstractGuidingFunction,Hilbert::AbstractHilbertSpace,Buffer = allocate_GWF_buffer(logψ,x))
    weights = zeros(length(H.moves))
    Hxx = get_diagonal(H)
    get_markov_weights!(weights,x,H,logψ,Hilbert,Buffer)
    # println(weights)
    return getLocalEnergy(x,weights,Hxx)
end

# function getLocalEnergy(WE::AbstractWalkerEnsemble,α,Hxx::DiagonalOperator)
#     Config = getConfig(WE,α)
#     moveWeights = getMoveWeights(WE,α)
#     return getLocalEnergy(Config,moveWeights,Hxx)
# end

function getLocalEnergyWalkers_before(Walkers::AbstractWalkerEnsemble,Hxx::DiagonalOperator)
    num = 0.
    denom = 0.
    WalkerWeights = getWalkerWeights(Walkers)
    localEnergies = getLocalEnergies(Walkers)

    num = WalkerWeights' * localEnergies
    denom = sum(WalkerWeights)
    
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
             Hilbert::AbstractHilbertSpace, parallelizer::AbstractParallelizationScheme, logger::AbstractLogger, RNG::Random.AbstractRNG)

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
- `logger::AbstractLogger`: The logger object used to record simulation progress.
- `RNG::Random.AbstractRNG`: The random number generator used for stochastic processes in the simulation.

# Description
This function performs the GFMC simulation by evolving the walker ensemble using the provided propagator and guiding function. It measures observables diagonal in the computational Basis at each step and applies reconfiguration to control the walker population. The simulation is parallelized according to the specified parallelization scheme.

# Notes
- The function modifies the `Walkers` object in place.
- Ensure that all input objects are properly initialized before calling this function.
"""
function runGFMC!(Walkers::AbstractWalkerEnsemble,Observables::AbstractObserver,reconfiguration::AbstractReconfigurationScheme,range,propagator::AbstractPropagator,logψ::AbstractGuidingFunction,H::AbstractSignFreeOperator,Hilbert::AbstractHilbertSpace,parallelizer::AbstractParallelizationScheme,logger::AbstractLogger,RNG::Random.AbstractRNG)
    compute_GWF_buffers!(Walkers,logψ)
    for i in range
        propagateWalkers!(Walkers,H,logψ,Hilbert,propagator,parallelizer,RNG)
        saveObservables_before!(Observables,i,Walkers,H,reconfiguration)
        reconfigurateWalkers!(Walkers,reconfiguration,RNG)
        saveObservables_after!(Observables,i,Walkers,H,reconfiguration)
        write_log(logger,i,range,Walkers,Observables,reconfiguration)
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
    Walkers::WE
    Propagator::Prop
    H::SFO
    logψ::GF
    Hilbert::HS
    parallelization::PS
    reconfiguration::RS
end

function GFMCProblem(config::AbstractConfig,NWalkers::Integer,prop::AbstractPropagator,H::AbstractSignFreeOperator,Hilbert::AbstractHilbertSpace,logψ::AbstractGuidingFunction;parallelization = MultiThreaded(num_tasks_default(NWalkers)),reconfiguration = MinimalReconfiguration(NWalkers))
    WE = allocate_walkerEnsemble(config,logψ,NWalkers,H)

    moves_vals = get_offdiagonal_elements(H)
    
    for (i) in eachindex(moves_vals)
        move = get_move(H,i)
        sites = affected_sites(move)
        in_bounds = checkbounds(Bool, config,sites)
        @assert in_bounds "move $i is out of bounds"
    end

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
- `logger`: Logger object to record simulation progress. Default is `NoLogger()`.
- `parallelization`: Parallelization strategy for the simulation, defaulting to `MultiThreaded(NWalkers)`, which tries to pick an optimal number of tasks to spawn.

# Returns
A `GFMCProblem` instance configured with the specified parameters.

# Notes
- The `logψ` and `logpsi` parameters are interchangeable and exactly one needs to be provided.
"""
function GFMCProblem(config::AbstractConfig,NWalkers::Integer,prop::AbstractPropagator; H, Hilbert, logψ = nothing, logpsi = nothing,kwargs...)
    
    if logψ !== nothing && logpsi !== nothing
        throw(ArgumentError("Only one of `logψ` or `logpsi` can be specified."))
    end

    GWF = logψ !== nothing ? logψ : logpsi
    return GFMCProblem(config,NWalkers,prop,H, Hilbert, GWF; kwargs...)
end

"""
    runGFMC!(prob::GFMCProblem, Observables::AbstractObserver, range; 
             logger = default_logger(), rng = Random.default_rng(), 
             reconfiguration = prob.reconfiguration, Propagator = prob.Propagator, 
             logψ = prob.logψ, H = prob.H, Hilbert = prob.Hilbert, 
             parallelization = prob.parallelization)

Run the Green's Function Monte Carlo (GFMC) simulation for the given problem.

# Arguments
- `prob::GFMCProblem`: The GFMC problem instance containing the system configuration and parameters.
- `Observables::AbstractObserver`: An observer object to track and record observables during the simulation.
- `range`: Integer or range: The range of iterations or steps over which the simulation will be performed.

# Keyword Arguments to override default values
- `logger`: (Optional) A logger instance for logging simulation progress. Defaults to `default_logger()`. Use `NoLogger()` or `nothing` to disable logging.
- `rng`: (Optional) A random number generator to ensure reproducibility. Defaults to `Random.default_rng()`.
- `reconfiguration`: (Optional) The reconfiguration method to be used during the simulation. Defaults to `prob.reconfiguration`.
- `Propagator`: (Optional) The propagator function for the simulation. Defaults to `prob.Propagator`.
- `logψ`: (Optional) The logarithm of the wavefunction. Defaults to `prob.logψ`.
- `H`: (Optional) The Hamiltonian operator. Defaults to `prob.H`.
- `Hilbert`: (Optional) The Hilbert space of the system. Defaults to `prob.Hilbert`.
- `parallelization`: (Optional) The parallelization strategy for the simulation. Defaults to `prob.parallelization`.

# Returns
This function modifies the `prob` and `Observables` in place to reflect the results of the simulation.

# Notes
- Ensure that the `prob` and `Observables` are properly initialized before calling this function.
- If `range` is provided as an integer, it will be converted to a range `1:range`.
- If `logger` is `nothing`, it will default to `NoLogger()`.
"""
function runGFMC!(prob::GFMCProblem,Observables::AbstractObserver,range; 
    logger = default_logger(),
    rng = Random.default_rng(),
    reconfiguration = prob.reconfiguration,
    Propagator = prob.Propagator,
    logψ = prob.logψ,
    H = prob.H,
    Hilbert = prob.Hilbert,
    parallelization = prob.parallelization)

    if range isa Integer
        range = 1:range
    end
    if isnothing(logger)
        logger = NoLogger()
    end
    runGFMC!(prob.Walkers,Observables,reconfiguration,range,Propagator,logψ,H,Hilbert,parallelization,logger,rng)
end

"""
    ProblemEnsemble{P<:AbstractGFMCProblem} <: AbstractGFMCProblem

Run multiple GFMC problems in parallel. Useful to estimate errors.

# Usage example:
```
P = ProblemEnsemble([GFMCProblem1, GFMCProblem2, ...])

problems = ProblemEnsemble([GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ, H, Hilbert) for _ in 1:10])

Observers = [ConfigObserver("output_\$i.h5",startConfig, NSteps, NWalkers) for i in 1:10] #note that each observer must have its own file

runGFMC!(problems, NoObserver(),100) #equilibrate
runGFMC!(problems, Observers)
```

"""
struct ProblemEnsemble{P<:AbstractGFMCProblem} <: AbstractGFMCProblem
    problems::Vector{P}
end

function runGFMC!(P::ProblemEnsemble,Observers,args...;kwargs...)
    Threads.@threads for i in eachindex(P.problems,Observers)
        prob = P.problems[i]
        Observer = Observers[i]
        runGFMC!(prob,Observer,args...;kwargs...)
    end
    return Observers
end

function runGFMC!(P::ProblemEnsemble,Observer::NoObserver,args...;kwargs...)
    runGFMC!(P::ProblemEnsemble,[NoObserver() for _ in P.problems],args...;kwargs...)
end

getConfigs(p::GFMCProblem) = getConfigs(p.Walkers)