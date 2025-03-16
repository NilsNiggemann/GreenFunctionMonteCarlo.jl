"""
    AbstractHilbertSpace

An abstract type representing a Hilbert space. This type serves as a base for defining various specific Hilbert spaces used in the project. Generally, a hilbert space should be defined by the number of sites, the number of local degrees of freedom (i.e. spin), and the constraint that the configurations must satisfy. 
# Interface
- `constraint(HilbertSpace::AbstractHilbertSpace)`: Return the constraint that the configurations in the Hilbert space must satisfy.
- `Base.size(HilbertSpace::AbstractHilbertSpace)`: Return the size of a config in the Hilbert space.
# Interface (optional)
- fulfills_constraint(x,HilbertSpace::AbstractHilbertSpace): Check if the given configuration `x` satisfies the constraint of the specified `HilbertSpace`. Defaults to `constraint(HilbertSpace)(x)`.
"""
abstract type AbstractHilbertSpace end

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