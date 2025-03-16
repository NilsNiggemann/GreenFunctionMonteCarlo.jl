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
    pre_move_affect!(Buffer,x,dx,logψ)
    Hxy = get_offdiagonal_elements(H)
    moves = get_moves(H)

    for i in eachindex(moves,weights)
        ψx´_ψx = exp(log_psi_diff(x,moves[i],logψ,Buffer,Hilbert))
        weights[i] = - ψx´_ψx * Hxy[i]
    end

    return weights
end

getLocalEnergy(weights::AbstractVector) = -sum(weights)
getLocalEnergy(x::AbstractConfig,weights::AbstractVector,Hxx::DiagonalOperator) = getLocalEnergy(weights) + Hxx(x)