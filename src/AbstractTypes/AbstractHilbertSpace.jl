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