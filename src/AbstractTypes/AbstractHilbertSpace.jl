"""
    AbstractHilbertSpace

An abstract type representing a Hilbert space. This type serves as a base for defining various specific Hilbert spaces used in the project. Generally, a hilbert space should be defined by the number of sites, the number of local degrees of freedom (i.e. Spin), and the constraints that the configurations must satisfy. 
"""
abstract type AbstractHilbertSpace end

"""
    AbstractConstraint

An abstract type representing a constraint in the context of a Hilbert space.
# Interface
- `(c)(x::AbstractArray)`: Check if the configuration `x` satisfies the constraint.
"""
abstract type AbstractConstraint end

constraints(::AbstractHilbertSpace) = error("No constraints defined for this Hilbert space.")
size(::AbstractHilbertSpace) = error("No size defined for this Hilbert space.")

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
    result = fulfills_constraint(x, HilbertSpace)
    apply_inverse!(x, move)
    return result
end

