"""
    AbstractHilbertSpace

An abstract type representing a Hilbert space. This type serves as a base for defining various specific Hilbert spaces used in the project. Generally, a hilbert space should be defined by the number of sites, the number of local degrees of freedom (i.e. Spin), and the constraints that the configurations must satisfy. 
# Interface
- `constraints(HilbertSpace::AbstractHilbertSpace)`: Return the constraints that the configurations in the Hilbert space must satisfy.
- `Base.size(HilbertSpace::AbstractHilbertSpace)`: Return the size of the Hilbert space.
- fulfills_constraint(x,HilbertSpace::AbstractHilbertSpace): Check if the given configuration `x` satisfies the constraints of the specified `HilbertSpace`.
"""
abstract type AbstractHilbertSpace end

"""
    AbstractConstraint

An abstract type representing a constraint in the context of a Hilbert space.
# Interface
- `(c)(x::AbstractArray)`: Check if the configuration `x` satisfies the constraint.
"""
abstract type AbstractConstraint end

"""get the constraints that the configurations in the Hilbert space must satisfy."""
function constraints end

"""
    fulfills_constraint(x::AbstractArray, HilbertSpace::AbstractHilbertSpace)

Check if the given configuration `x` satisfies the constraints of the specified `HilbertSpace`.

# Arguments
- `x::AbstractArray`: The configuration to be checked.
- `HilbertSpace::AbstractHilbertSpace`: The Hilbert space whose constraints need to be satisfied.

# Returns
- `Bool`: `true` if the configuration satisfies the constraints, `false` otherwise.
"""
function fulfills_constraint end
