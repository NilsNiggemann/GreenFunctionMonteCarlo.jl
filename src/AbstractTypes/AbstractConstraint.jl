"""
    AbstractConstraint

An abstract type representing a constraint in the context of a Hilbert space.
# Interface
- `(c)(x::AbstractArray)`: Check if the configuration `x` satisfies the constraint.
"""
abstract type AbstractConstraint end

struct MultipleConstraints{Cons<:Tuple} <: AbstractConstraint
    constraints::Cons
end
(c::MultipleConstraints)(x) = all(c_i(x) for c_i in c.constraints)

struct NoConstraint <: AbstractConstraint end
(cons::NoConstraint)(x) = true

"""
    fulfills_constraint(x::AbstractArray, HilbertSpace::AbstractHilbertSpace)

Check if the given configuration `x` satisfies the constraint of the specified `HilbertSpace`.

# Arguments
- `x::AbstractArray`: The configuration to be checked.
- `HilbertSpace::AbstractHilbertSpace`: The Hilbert space whose constraint need to be satisfied.

# Returns
- `Bool`: `true` if the configuration satisfies the constraint, `false` otherwise.
"""
function fulfills_constraint(x::AbstractConfig,H::AbstractHilbertSpace) 
    c = constraint(H)
    return c(x)
end