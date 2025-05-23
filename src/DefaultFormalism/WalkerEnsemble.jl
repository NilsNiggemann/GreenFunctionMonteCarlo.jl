struct WalkerEnsemble{ConfType<:AbstractConfig,GWB<:AbstractGuidingFunctionBuffer} <: AbstractWalkerEnsemble
    Configs::Vector{ConfType}
    WalkerWeights::Vector{Float64}
    MoveWeights::Vector{Vector{Float64}}
    Buffers::Vector{GWB}
    reconfigurationList::Vector{Int}
end

getConfig(X::WalkerEnsemble, α) = X.Configs[α]
getConfigs(X::WalkerEnsemble) = X.Configs
getMoveWeights(X::WalkerEnsemble, α) = X.MoveWeights[α]
getWalkerWeights(X::WalkerEnsemble) = X.WalkerWeights
getBuffer(X::WalkerEnsemble, α) = X.Buffers[α]
getReconfigurationList(X::WalkerEnsemble) = X.reconfigurationList

function allocate_walkerEnsemble(conf, logψ::AbstractGuidingFunction,NWalkers::Integer,numMoves::Integer)
    configs = [copy(conf) for _ in 1:NWalkers]
    weights = zeros(NWalkers)
    move_weights = [zeros(numMoves) for _ in 1:NWalkers]
    buffers = allocate_GWF_buffers(logψ, conf, NWalkers)
    reconfigurationList = zeros(Int, NWalkers)
    return WalkerEnsemble(configs, weights, move_weights, buffers, reconfigurationList)
end
allocate_walkerEnsemble(conf, logψ::AbstractGuidingFunction,NWalkers::Integer,H::AbstractSignFreeOperator) = allocate_walkerEnsemble(conf, logψ, NWalkers, length(get_offdiagonal_elements(H)))