struct zeroDiagOperator <: AbstractOperator end

struct LocalOperator{MoveType,T,DiagType} <: AbstractSignFreeOperator
    moves::SparseArrays.SparseMatrixCSC{MoveType,Int}
    weights::Vector{T}
    diag::DiagType
end

struct SparseMove{MoveType} <: AbstractMove
    dx::MoveType
end

get_diagonal(O::LocalOperator) = O.diag
get_offdiagonal_elements(O::LocalOperator) = O.weights
get_move(O::LocalOperator,idx::Integer) = SparseMove(@view O.moves[:,idx])
apply!(x::AbstractConfig, move::SparseMove) = x .+= move