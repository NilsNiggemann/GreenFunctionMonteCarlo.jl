struct ZeroDiagOperator <: AbstractOperator end
ZeroDiagOperator(x::AbstractConfig) = 0.

struct SparseMove{T} <: AbstractMove
    dx::SparseArrays.SparseVector{T,Int}
end

@inline function apply!(x::AbstractConfig, move::SparseMove)
    for i in move.dx.nzind
        x[i] += move.dx[i]
    end
end
@inline apply_inverse!(x::AbstractConfig, move::SparseMove) = x .-= move.dx
@inline function apply!(x::AbstractConfig{Bool}, move::SparseMove{Bool})
    for i in move.dx.nzind
        x[i] = !x[i]
    end
end
@inline apply_inverse!(x::AbstractConfig{Bool}, move::SparseMove{Bool}) = apply!(x, move)
struct LocalOperator{MoveType,T,DiagType} <: AbstractSignFreeOperator
    moves::Vector{MoveType}
    weights::Vector{T}
    diag::DiagType
end
@inline get_move(O::LocalOperator{<:Any},idx::Integer) = O.moves[idx]

isapplicable(x::AbstractConfig{Bool}, move::SparseMove{Bool}, c::HardCoreConstraint) = true

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
