struct WalkerEnsemble{ConfType<:AbstractConfig,GWB<:AbstractGuidingFunctionBuffer} <: AbstractWalkerEnsemble
    Configs::Vector{ConfType}
    WalkerWeights::Vector{Float64}
    MoveWeights::Vector{Vector{Float64}}
    Buffers::Vector{GWB}
end

getConfig(X::WalkerEnsemble, α) = X.Configs[α]
getMoveWeights(X::WalkerEnsemble, α) = X.MoveWeights[α]
getWalkerWeights(X::WalkerEnsemble) = X.WalkerWeights
getBuffer(X::WalkerEnsemble, α) = X.Buffers[α]

function allocate_walkerEnsemble(conf::ConfType, logψ::AbstractGuidingFunction,NWalkers::Integer,numMoves::Integer) where {ConfType<:AbstractConfig}
    configs = [copy(conf) for _ in 1:NWalkers]
    weights = zeros(NWalkers)
    move_weights = [zeros(numMoves) for _ in 1:NWalkers]
    buffers = allocate_GWF_buffers(logψ, conf, NWalkers)
    return WalkerEnsemble(configs, weights, move_weights, buffers)
end