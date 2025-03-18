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
    for α in eachindex(WalkerWeights,Walkers)
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

struct NoObserver <: AbstractObserver end
saveObservables_before!(::NoObserver,i,Walkers,propagator,reconfigurations) = nothing
saveObservables_after!(::NoObserver,i,Walkers,propagator,reconfigurations) = nothing

function runGFMC!(Walkers::AbstractWalkerEnsemble,Observables::AbstractObserver,reconfiguration::AbstractReconfigurationScheme,range,propagator::AbstractPropagator,logψ::AbstractGuidingFunction,H::AbstractSignFreeOperator,Hilbert::AbstractHilbertSpace,parallelizer::AbstractParallelizationScheme,RNG::Random.AbstractRNG)
    iter = 0
    for i in range
        iter += 1
        propagateWalkers!(Walkers,H,logψ,Hilbert,propagator,parallelizer,RNG)
        saveObservables_before!(Observables,i,Walkers,propagator,reconfiguration)
        reconfigurateWalkers!(Walkers,reconfiguration,RNG)
        saveObservables_after!(Observables,i,Walkers,propagator,reconfiguration)

        if iter%1000 == 0 # recompute buffers only occasionally to avoid accumulation of floating point errors 
            fill_all_Buffers!(prob,nThreads)
        end
    end
    return Observables
end

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

function GFMCProblem(config::AbstractConfig,NWalkers::Integer,prop::AbstractPropagator; H, Hilbert, logψ = nothing, logpsi = nothing,parallelization = MultiThreaded(NWalkers))
    
    if logψ !== nothing && logpsi !== nothing
        throw(ArgumentError("Only one of `logψ` or `logpsi` can be specified."))
    end

    GWF = logψ !== nothing ? logψ : logpsi
    return GFMCProblem(config,NWalkers,prop,H, Hilbert, GWF; parallelization)
end

runGFMC!(prob::GFMCProblem,Observables::AbstractObserver,range, rng = Random.default_rng()) = runGFMC!(prob.WE,Observables,prob.reconfiguration,range,prob.Propagator,prob.logψ,prob.H,prob.Hilbert,prob.parallelization,rng)
runGFMC!(prob::GFMCProblem,Observables::AbstractObserver,NSteps::Integer, rng = Random.default_rng()) = runGFMC!(prob,Observables,1:NSteps,rng)