struct ZeroDiagOperator <: AbstractOperator end

struct SparseMove{T} <: AbstractMove
    dx::SparseArrays.SparseVector{T,Int}
end


@inline apply!(x::AbstractConfig, move::SparseMove) = x .+= move.dx
@inline apply_inverse!(x::AbstractConfig, move::SparseMove) = x .-= move.dx
@inline apply!(x::AbstractConfig{Bool}, move::SparseMove{Bool}) = x .= x .⊻ move.dx
@inline apply_inverse!(x::AbstractConfig{Bool}, move::SparseMove{Bool}) = apply!(x, move)
struct LocalOperator{MoveType,T,DiagType} <: AbstractSignFreeOperator
    moves::Vector{MoveType}
    weights::Vector{T}
    diag::DiagType
end
@inline get_move(O::LocalOperator{<:Any},idx::Integer) = O.moves[idx]

_move_type(::LocalOperator{MoveType}) where {MoveType} = MoveType

function localOperator(moves::AbstractVector, weights, diag)
    @assert all(>=(0), weights) "Weights must be strictly nonnegative"
    @assert length(moves) == length(weights) "Number of columns of moves must match number of weights"
    
    any(iszero, weights) && @warn "Some weights are zero"
    
    spmoves = [SparseMove(SparseArrays.sparse(dx)) for dx in moves]
    return LocalOperator(spmoves, weights, diag)
end
function localOperator(moves::AbstractMatrix, weights, diag)
    localOperator(eachcol(moves), weights, diag)
end
get_diagonal(O::LocalOperator) = O.diag
get_offdiagonal_elements(O::LocalOperator) = O.weights
