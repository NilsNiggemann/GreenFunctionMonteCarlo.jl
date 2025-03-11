"""
    AbstractWalkerEnsemble{T,N}

An abstract type that represents an ensemble of configurations (i.e. spins, bosons on a lattice for each walker) in the context of the Green Function Monte Carlo project. 
This type is a subtype of `AbstractArray{T,N}`, indicating that it behaves like an N-dimensional array with elements of type `T`.

# Parameters
- `T`: The type of the elements in the configurations.
- `N`: The dimensionality given by the dimension of the system plus an additional dimension for the walkers.

# Interface
- `parent(X::AbstractWalkerEnsemble)`: return the parent array
- apply_move!(X::AbstractWalkerEnsemble, move::Any): apply a move to the configuration
- fulfills_constraints(X::AbstractWalkerEnsemble, HilbertSpace::AbstractHilbertSpace): check if the configuration satisfies the constraints of the Hilbert space
- getConfig(X::AbstractWalkerEnsemble,WalkerIndex): get the configuration at index `WalkerIndex` as a view.
- eachConfig(X::AbstractWalkerEnsemble): iterate over the configurations of the ensemble.
"""
abstract type AbstractWalkerEnsemble{T,N} <: AbstractArray{T,N} end

parent(X::AbstractWalkerEnsemble) = throw(MethodError(parent, (X,)))

size(X::AbstractWalkerEnsemble) = size(parent(X))
getindex(X::AbstractWalkerEnsemble, i::Int) = getindex(parent(X), i)
getindex(X::AbstractWalkerEnsemble, I::Vararg{Int,N}) = getindex(parent(X), I...)
getindex(X, I...) = getindex(parent(X), I...)

IndexStyle(::Type{<:AbstractWalkerEnsemble}) = IndexStyle(parenttype(AbstractWalkerEnsemble))

setindex!(X::AbstractWalkerEnsemble, v, i::Int) = setindex!(parent(X), v, i)
setindex!(X::AbstractWalkerEnsemble, v, I::Vararg{Int,N}) = setindex!(parent(X), v, I...)
setindex!(X, I...) = setindex!(parent(X), X, I...)

getConfig(X::AbstractWalkerEnsemble{T,N},WalkerIndex) where {T,N} = selectdim(X,N,WalkerIndex)
eachConfig(X::AbstractWalkerEnsemble{T,N}) where {T,N} = eachslice(X,dims=N)

"""
    fulfills_constraint(x::AbstractArray, HilbertSpace::AbstractHilbertSpace)

Check if the given configuration `x` satisfies the constraints of the specified `HilbertSpace`.

# Arguments
- `x::AbstractArray`: The configuration to be checked.
- `HilbertSpace::AbstractHilbertSpace`: The Hilbert space whose constraints need to be satisfied.

# Returns
- `Bool`: `true` if the configuration satisfies the constraints, `false` otherwise.
"""
function fulfills_constraint(x::AbstractArray,HilbertSpace::AbstractHilbertSpace)
    constraints = constraints(HilbertSpace)
    return all(c(X) for c in constraints)
end

"""
    struct InverseMove{T}

A type representing an inverse move in a Monte Carlo simulation. The `InverseMove` struct is used to encapsulate the concept of an inverse operation corresponding to a given move in the simulation.
"""
struct InverseMove{T}
    x::T
end

"""
    apply_move!(x::AbstractArray, move::AbstractArray)

Applies a move to the configuration `x`.

# Arguments
- `x::AbstractArray`: An instance of a type that is a subtype of `AbstractArray`.
- `move::AbstractArray`: The move to be applied to the configuration.
"""
apply_move!(x::AbstractArray, move::Any) = throw(MethodError(apply_move!, (x, move)))
apply_move!(x::AbstractArray, move::AbstractArray) = x .+= move
apply_move!(x::AbstractArray, invmove::InverseMove) = apply_inverse!(x, invmove.x)
apply_inverse!(x::AbstractArray, move::AbstractArray) = x .-= move

"""
    isapplicable(x::AbstractArray, move::AbstractArray, HilbertSpace::AbstractHilbertSpace)

Check if a given move is applicable to the current configuration within a specified Hilbert space.

# Arguments
- `x::AbstractArray`: The current configuration.
- `move::AbstractArray`: The proposed move to be applied.
- `HilbertSpace::AbstractHilbertSpace`: The Hilbert space in which the configuration and move are defined.

# Returns
- `Bool`: `true` if the move is applicable, `false` otherwise.

Defaults to applying the move to the configuration, checking if the constraints are satisfied, and then reverting the move.

Make sure to implement this function for specific subtypes of `AbstractArray` and `AbstractHilbertSpace` to ensure optimal performance.
"""
function isapplicable(x::AbstractArray, move::AbstractArray, HilbertSpace::AbstractHilbertSpace)
    apply_move!(x, move)
    result = fulfills_constraints(x, HilbertSpace)
    apply_inverse!(x, move)
    return result
end

