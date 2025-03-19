struct ZeroDiagOperator <: DiagonalOperator end
(::ZeroDiagOperator)(x::AbstractConfig) = 0.

struct SparseMove{T,V1<:AbstractVector{Int},V2<:AbstractVector{T}} <: AbstractMove
    inds::V1
    vals::V2
end

@inline function apply!(x::AbstractConfig, move::SparseMove)
    for (i,dx) in zip(move.inds, move.vals)
        x[i] += dx
    end
end

@inline function apply_inverse!(x::AbstractConfig, move::SparseMove)
    for (i,dx) in zip(move.inds, move.vals)
        x[i] -= dx
    end
end

struct FlipMove{Vec} <: AbstractMove
    inds::Vec
end

@inline function apply!(x::AbstractConfig{Bool}, move::FlipMove)
    for i in move.inds
        x[i] = !x[i]
    end
end
@inline apply_inverse!(x::AbstractConfig{Bool}, move::FlipMove) = apply!(x, move)

struct LocalOperator{MoveType,T,DiagType} <: AbstractSignFreeOperator
    moves::Vector{MoveType}
    off_diag::Vector{T}
    diag::DiagType
end

@inline get_move(O::LocalOperator{<:Any}, idx::Integer) = O.moves[idx]
@inline get_diagonal(O::LocalOperator) = O.diag
@inline get_offdiagonal_elements(O::LocalOperator) = O.off_diag

isapplicable(x::AbstractConfig{Bool}, move::FlipMove, c::HardCoreConstraint) = true

_move_type(::LocalOperator{MoveType}) where {MoveType} = MoveType

function localOperator(moves::AbstractVector, off_diag::AbstractVector, diag, H::AbstractHilbertSpace)
    @assert length(moves) == length(off_diag) "Number of moves must match number of Off-diagonal elements"
    @assert all(<=(0), off_diag) "Off-diagonal elements must be strictly negative"
    numsites = size(H)
    @assert all(length.(moves) .== numsites) "sizes of moves ($(length.(moves))) must correspond to number of sites $numsites"
    
    any(iszero, off_diag) && @warn "Some weights are zero"
    
    num_affected = [sum(≠(0), m) for m in moves]

    if all(==(num_affected[1]), num_affected)
        numsites = Val(num_affected[1])
        spmoves = [_convert_to_sparse_move_SVec(move,numsites) for move in moves]
    else
        max_num = Val(maximum(num_affected))
        spmoves = [_convert_to_sparse_move_SmallVec(move,max_num) for move in moves]
    end
    return LocalOperator(spmoves, off_diag, diag)
end

function _convert_to_sparse_move_SVec(dx::AbstractVector{T},::Val{Nsites}) where {T,Nsites}
    nzinds = findall(!iszero,dx)

    inds = SA.SVector{Nsites,Int}(nzinds)
    vals = SA.SVector{Nsites,T}(dx[nzinds])
    move = SparseMove(inds,vals)
    return _maybe_reduce_to_flipMove(move)
end
function _convert_to_sparse_move_SmallVec(dx::AbstractVector{T},::Val{Nsites}) where {T,Nsites}
    nzinds = findall(!iszero,dx)

    inds = SmallCollections.SmallVector{Nsites,Int}(nzinds)
    vals = SmallCollections.SmallVector{Nsites,T}(dx[nzinds])
    move = SparseMove(inds,vals)
    return _maybe_reduce_to_flipMove(move)
end

function _maybe_reduce_to_flipMove(move::SparseMove{Bool,Nsites}) where {Nsites}
    return FlipMove(move.inds)
end
_maybe_reduce_to_flipMove(move::Any) = move