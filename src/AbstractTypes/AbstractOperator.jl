"""
    AbstractOperator

An abstract type representing a general operator. This serves as a base type for defining various specific operators in the context of the Green Function Monte Carlo project.
"""
abstract type AbstractOperator end
abstract type OffdiagonalOperator end

"""
    DiagonalOperator

An abstract type representing a diagonal operator in the context of Green Function Monte Carlo simulations.
A diagonal operator is special in the sense that it will not change the configuration of the system when applied to it and will only return a number.
# Interface: 
- (D::DiagonalOperator)(x) : return the value of the operator applied to the configuration `x`
- (D::DiagonalOperator)(x, params): return the value of the operator applied to the configuration `x` with parameters `params`, which can be used to store buffers.
"""
abstract type DiagonalOperator end

"""
    get_move(O::AbstractOperator,idx::Integer)
Return the move associated with the operator `O` at index `idx`. 
"""
function get_move end

function (D::DiagonalOperator) end

"""
    AbstractSignFreeOperator <: AbstractOperator

An abstract type representing a sign-free operator in the context of Green Function Monte Carlo simulations. 
# Interface:
- get_diagonal(O::AbstractSignFreeOperator): return the diagonal operator associated with the sign-free operator `O`
- get_offdiagonal_elements(O::AbstractSignFreeOperator): return the weights associated with the off-diagonal operator `O`
"""
abstract type AbstractSignFreeOperator <: AbstractOperator end
"""
    get_offdiagonal_elements(O::AbstractSignFreeOperator)

Return the weights associated with an `OffdiagonalOperator` object `O`. 
"""
function get_offdiagonal_elements end

function get_diagonal end

abstract type AbstractMove end

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
function apply! end
apply!(C::AbstractConfig, move::InverseMove) = apply_inverse!(C, move.x)

"""
    isapplicable(x::AbstractConfig, move::AbstractArray, HilbertSpace::AbstractHilbertSpace)

Check if a given move is applicable to the current configuration within a specified Hilbert space.

# Arguments
- `x::AbstractConfig`: The current configuration.
- `move::AbstractMove`: The proposed move to be applied.
- `HilbertSpace::AbstractHilbertSpace`: The Hilbert space in which the configuration and move are defined.

# Returns
- `Bool`: `true` if the move is applicable, `false` otherwise.

Defaults to applying the move to the configuration, checking if the constraints are satisfied, and then reverting the move.

Make sure to implement this function for specific subtypes of `AbstractConfig` and `AbstractHilbertSpace` to ensure optimal performance.
"""
function isapplicable(x::AbstractConfig, move::AbstractMove, HilbertSpace::AbstractHilbertSpace)
    apply!(x, move)
    result = fulfills_constraints(x, HilbertSpace)
    apply_inverse!(x, move)
    return result
end

