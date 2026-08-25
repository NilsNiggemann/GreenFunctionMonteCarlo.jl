"""
    LocalOffdiagonalObservable{MoveType,WeightFunc} <: AbstractOffdiagonalObservable

A generic, closure-based `AbstractOffdiagonalObservable`: a fixed list of candidate moves plus a
user-supplied `weightfunc(idx, ψratio, x) -> Real` computing the non-negative sampling weight for
each move. Construct via [`offdiagonalObservable`](@ref).
"""
struct LocalOffdiagonalObservable{MoveType,WeightFunc} <: AbstractOffdiagonalObservable
    moves::Vector{MoveType}
    weightfunc::WeightFunc
end

@inline get_move(O::LocalOffdiagonalObservable, idx::Integer) = O.moves[idx]
@inline n_moves(O::LocalOffdiagonalObservable) = length(O.moves)
@inline observable_weight(O::LocalOffdiagonalObservable, idx::Integer, ψratio::Real, x::AbstractConfig) = O.weightfunc(idx, ψratio, x)

"""
    offdiagonalObservable(moves::AbstractVector, weightfunc, Hilbert::AbstractHilbertSpace)

Construct a [`LocalOffdiagonalObservable`](@ref) from a list of candidate `moves` (in the same
dense-vector-per-move format accepted by `localOperator`) and a per-move sampling-weight
function `weightfunc(idx, ψratio, x) -> Real`, where `ψratio = ψ(x′)/ψ(x)` is the guiding-function
ratio for that move.

Unlike `localOperator`, this does **not** assert any sign constraint on the resulting
weights: `weightfunc` is responsible for returning a non-negative value for every move (this is
required for the categorical sampling in `sample_and_apply_observable!` to be well-defined). A
common way to guarantee this for a Hermitian off-diagonal operator is a `cos²` phase weighting
(e.g. for a structure factor `B(q)`), which avoids a sign problem without constraining the
underlying physical matrix elements.
"""
function offdiagonalObservable(moves::AbstractVector, weightfunc, Hilbert::AbstractHilbertSpace)
    numsites = size(Hilbert)
    @assert all(length.(moves) .== numsites) "sizes of moves ($(length.(moves))) must correspond to number of sites $numsites"

    num_affected = [sum(≠(0), m) for m in moves]

    spmoves = if all(==(num_affected[1]), num_affected)
        Nsites = Val(num_affected[1])
        [_convert_to_sparse_move_SVec(move, Nsites) for move in moves]
    else
        max_num = Val(maximum(num_affected))
        [_convert_to_sparse_move_SmallVec(move, max_num) for move in moves]
    end
    return LocalOffdiagonalObservable(spmoves, weightfunc)
end

# Mirrors `_get_markov_weights!` (ManyWalkerGFMC.jl) structurally: same pre_move_affect! + per-move
# log_psi_diff loop, but calls the operator's own `observable_weight` instead of the fixed
# `-Hxy*exp(Δ)` propagation formula.
function get_observable_weights!(weights::AbstractVector, x::AbstractConfig, O::AbstractOffdiagonalObservable,
                                  logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, Buffer)
    pre_move_affect!(Buffer, x, logψ)
    for i in eachindex(weights)
        move = get_move(O, i)
        Δ = log_psi_diff(x, move, logψ, Buffer, Hilbert)
        ψratio = exp(Δ)
        weights[i] = observable_weight(O, i, ψratio, x)
    end
    return weights
end

# Mirrors `performMarkovStep!` (ManyWalkerGFMC.jl) structurally, additionally returning `sum(weights)`
# as the physical seed weight (SpiderWebModel's `wsum`/`Operator_weight`). Permanently mutates `x`.
function sample_and_apply_observable!(x::AbstractConfig, weights::AbstractVector, O::AbstractOffdiagonalObservable, rng::Random.AbstractRNG)
    w = sum(weights)
    iszero(w) && return nothing, w
    moveidx = StatsBase.sample(rng, StatsBase.Weights(weights, w))
    move = get_move(O, moveidx)
    apply!(x, move)
    return move, w
end
