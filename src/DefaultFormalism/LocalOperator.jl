struct zeroDiagOperator <: AbstractOperator end

struct LocalOperator{MoveType,T,DiagType} <: AbstractSignFreeOperator
    moves::SparseArrays.SparseMatrixCSC{MoveType,Int}
    weights::Vector{T}
    diag::DiagType
end

get_diagonal(O::LocalOperator) = O.diag
get_offdiagonal_elements(O::LocalOperator) = O.weights
get_moves(O::LocalOperator) = O.moves