function get_markov_weights!(weights::AbstractVector,x::AbstractConfig,H::AbstractSignFreeOperator,logψ::AbstractGuidingFunction,Hilbert::AbstractHilbertSpace,Buffer::AbstractGuidingFunctionBuffer)
    newBuff = pre_move_affect!(Buffer,x,dx,logψ)
    Hxy = get_offdiagonal_elements(H)
    moves = get_moves(H)

    for i in eachindex(moves,weights)
        ψx´_ψx = exp(log_psi_diff(x,moves[i],logψ,newBuff,Hilbert))
        weights[i] = - ψx´_ψx * Hxy[i]
    end

    return weights
end