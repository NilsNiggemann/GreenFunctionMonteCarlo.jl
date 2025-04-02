"""
    ZeroDiagOperator <: DiagonalOperator

A struct representing a diagonal operator `H_xx = 0`
"""
struct ZeroDiagOperator <: DiagonalOperator end
(::ZeroDiagOperator)(x::AbstractConfig) = 0.

"""
    DiagOperator{F} <: DiagonalOperator

A type that can hold an arbitrary function `H_xx = f(x)` as a diagonal operator.
"""
struct DiagOperator{F} <: DiagonalOperator
    f::F
end
(F::DiagOperator)(x::AbstractConfig) = F.f(x)

"""
    OneBodyDiagOperator{T<:AbstractVector} <: DiagonalOperator

A type that represents arbitrary diagonal one-body interactions `H_xx = sum_{i} m_i x_i`
"""
struct OneBodyDiagOperator{T<:AbstractVector} <: DiagonalOperator
    m_i::T
end
(Hxx::OneBodyDiagOperator)(x::AbstractConfig) =  LinearAlgebra.dot(Hxx.m_i, x)

"""
    TwoBodyDiagOperator{T<:AbstractMatrix} <: DiagonalOperator

A type that represents arbitrary diagonal two-body interactions `H_xx = sum_{i,j} v_ij x_i x_j`
"""
struct TwoBodyDiagOperator{T<:AbstractMatrix} <: DiagonalOperator
    v_ij::T
end
(Hxx::TwoBodyDiagOperator)(x::AbstractConfig) =  LinearAlgebra.dot(x,Hxx.v_ij, x)


struct DiagOperatorSum{T<:Tuple} <: DiagonalOperator
    terms::T
end
(O::DiagOperatorSum)(x::AbstractConfig) = sum(term(x) for term in O.terms)

Base.:+(A::DiagonalOperator, B::DiagonalOperator) = DiagOperatorSum((A,B))
Base.:+(A::DiagOperatorSum, B::DiagonalOperator) = DiagOperatorSum((A...,B))
Base.:+(A::DiagonalOperator, B::DiagOperatorSum) = DiagOperatorSum((A,B...))
Base.:+(A::DiagOperatorSum, B::DiagOperatorSum) = DiagOperatorSum((A...,B...))

struct SparseMove{T,V1<:AbstractVector{Int},V2<:AbstractVector{T}} <: AbstractMove
    inds::V1
    vals::V2
end
@inline affected_sites(move::SparseMove) = move.inds
@inline move_dx(move::SparseMove) = move.vals

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

@inline affected_sites(move::FlipMove) = move.inds

@inline function move_dx_before(move::FlipMove,x::AbstractConfig{Bool})
    map(move.inds) do i
        !x[i] - x[i]
    end
end
@inline move_dx_after(move::FlipMove,x::AbstractConfig{Bool}) = -move_dx_before(move,x)

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
isapplicable(x::AbstractConfig{Bool}, move::FlipMove, c::OccupationNumberConstraint) = true

function isapplicable(x::AbstractConfig{<:Integer}, move::SparseMove, c::OccupationNumberConstraint)
    for (i,site) in enumerate(move.inds)
        if !(c.min_occupation <= x[site] + move.vals[i] <= c.max_occupation)
            return false
        end
    end
    return true
end

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