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

function runGFMC!(WE::AbstractWalkerEnsemble,Observables::AbstractObservables,propagator::AbstractPropagator,logψ::AbstractGuidingFunction,H::AbstractSignFreeOperator,Hilbert::AbstractHilbertSpace,parallelizer::AbstractParallelizationScheme,RNG::Random.AbstractRNG = Random.default_rng())
    reconfigurationList = getReconfigurationList(WE)
    iter = 0
    for i in range
        iter += 1
        propagateWalkers!(WE,H,logψ,Hilbert,propagator,w_avg_estimate,parallelizer,RNG)
        updateEnergies!(Observables,i,Walkers,weights,method)
        if reconfigure
            reconfiguration!(Walkers,Guiding_function_buffer,reconfigurationList,reconfiguration_buffer,weights)
        end
        saveObservables!(Observables,i,Walkers)

        if iter%1000 == 0 # recompute buffers only occasionally to avoid accumulation of floating point errors 
            fill_all_Buffers!(prob,nThreads)
        end
    end
    return Observables
end