function get_markov_weights!(weights::AbstractVector,x::AbstractConfig,H::AbstractSignFreeOperator,logψ::AbstractGuidingFunction,Hilbert::AbstractHilbertSpace,Buffer::AbstractGuidingFunctionBuffer)
    newBuff = pre_move_affect!(Buffer,x,dx,logψ)
    Hxy = get_offdiagonal_elements(H)
    moves = get_moves(H)

    for move in moves
        weight = log_psi_diff(x,move,logψ,newBuff,Hilbert)
        push!(weights,weight)
    end

    # LoopVectorization.@turbo for i in eachindex(weights)
    for i in eachindex(weights)
        weights[i] = exp(weights[i]) * Hxy[i]
    end
    return weights
end