"""
    AbstractConfig{T,N}

An abstract type that represents a configuration (i.e. spins, bosons on a lattice) in the context of the Green Function Monte Carlo project. 
This type is a subtype of `AbstractArray{T,N}`, indicating that it behaves like an N-dimensional array with elements of type `T`.

# Parameters
- `T`: The type of the elements in the configuration.
- `N`: The number of dimensions of the configuration.

# Interface
- `parent(x::AbstractConfig)`: return the parent array
- apply!(x::AbstractConfig, move::Any): apply a move to the configuration
- fulfills_constraints(x::AbstractConfig, HilbertSpace::AbstractHilbertSpace): check if the configuration satisfies the constraints of the Hilbert space
"""
abstract type AbstractConfig{T,N} <: AbstractArray{T,N} end

"""
    AbstractMove

An abstract type representing a move in the context of a Monte Carlo simulation. 
Subtypes of `AbstractMove` should implement specific types of moves or updates 
that can be applied to the configuration of the system being simulated.
"""
abstract type AbstractMove end


parent(x::AbstractConfig) = throw(MethodError(parent, (x,)))

size(A::AbstractConfig) = size(parent(A))
getindex(A::AbstractConfig, i::Int) = getindex(parent(A), i)
getindex(A::AbstractConfig, I::Vararg{Int,N}) = getindex(parent(A), I...)
getindex(A, I...) = getindex(parent(A), I...)

IndexStyle(::Type{<:AbstractConfig}) = IndexStyle(parenttype(AbstractConfig))

setindex!(A::AbstractConfig, v, i::Int) = setindex!(parent(A), v, i)
setindex!(A::AbstractConfig, v, I::Vararg{Int,N}) = setindex!(parent(A), v, I...)
setindex!(A, X, I...) = setindex!(parent(A), X, I...)

"""
    fulfills_constraints(x::AbstractConfig, HilbertSpace::AbstractHilbertSpace)

Check if the given configuration `x` satisfies the constraints of the specified `HilbertSpace`.

# Arguments
- `x::AbstractConfig`: The configuration to be checked.
- `HilbertSpace::AbstractHilbertSpace`: The Hilbert space whose constraints need to be satisfied.

# Returns
- `Bool`: `true` if the configuration satisfies the constraints, `false` otherwise.
"""
function fulfills_constraints(x::AbstractConfig,HilbertSpace::AbstractHilbertSpace)
    constraints = constraints(HilbertSpace)
    return all(c(x) for c in constraints)
end

"""
    struct InverseMove{T<:AbstractMove} <: AbstractMove

A type representing an inverse move in a Monte Carlo simulation. This type is parameterized by `T`, which must be a subtype of `AbstractMove`. The `InverseMove` struct is used to encapsulate the concept of an inverse operation corresponding to a given move in the simulation.

# Parameters
- `T`: A subtype of `AbstractMove` that specifies the type of move for which this struct represents the inverse.
"""
struct InverseMove{T<:AbstractMove} <: AbstractMove
    x::T
end
Base.inv(move::AbstractMove) = InverseMove(move)

"""
    apply!(x::AbstractConfig, move::AbstractMove)

Throws a `MethodError` indicating that the `apply!` function has not been implemented for the given `AbstractConfig` type and `move`.

# Arguments
- `x::AbstractConfig`: An instance of a type that is a subtype of `AbstractConfig`.
- `move::AbstractMove`: A move or operation to be applied to the configuration `x`.

# Throws
- `MethodError`: Always thrown to indicate that the method needs to be implemented for specific subtypes of `AbstractConfig`.
"""
apply!(x::AbstractConfig, move) = throw(MethodError(apply!, (x, move)))
apply!(C::AbstractConfig, move::AbstractArray) = C .+= move
apply!(C::AbstractConfig, move::InverseMove) = apply_inverse!(C, move.x)
apply_inverse!(C::AbstractConfig, move::AbstractArray) = C .-= move

"""
    isapplicable(x::AbstractConfig, move::AbstractArray, HilbertSpace::AbstractHilbertSpace)

Check if a given move is applicable to the current configuration within a specified Hilbert space.

# Arguments
- `x::AbstractConfig`: The current configuration.
- `move::AbstractArray`: The proposed move to be applied.
- `HilbertSpace::AbstractHilbertSpace`: The Hilbert space in which the configuration and move are defined.

# Returns
- `Bool`: `true` if the move is applicable, `false` otherwise.

Defaults to applying the move to the configuration, checking if the constraints are satisfied, and then reverting the move.

Make sure to implement this function for specific subtypes of `AbstractConfig` and `AbstractHilbertSpace` to ensure optimal performance.
"""
function isapplicable(x::AbstractConfig, move::AbstractArray, HilbertSpace::AbstractHilbertSpace)
    apply!(x, move)
    result = fulfills_constraints(x, HilbertSpace)
    apply_inverse!(x, move)
    return result
end

