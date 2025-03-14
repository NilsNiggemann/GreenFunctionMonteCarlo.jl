"""
    AbstractConfig{T,N}

An abstract type that represents a configuration (i.e. spins, bosons on a lattice). 
This type is a subtype of `AbstractArray{T,N}`.
# Parameters
- `T`: The type of the elements in the configuration.
- `N`: The number of dimensions of the configuration.

# Interface
- `Base.parent(x::AbstractConfig)`: return the parent array
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

Base.size(A::AbstractConfig) = size(parent(A))
Base.getindex(A::AbstractConfig, i::Int) = getindex(parent(A), i)
Base.getindex(A::AbstractConfig, I::Vararg{Int,N}) where N = getindex(parent(A), I...)
Base.getindex(A, I...) = getindex(parent(A), I...)

Base.IndexStyle(A::AbstractConfig) = IndexStyle(parent(A))

Base.setindex!(A::AbstractConfig, v, i::Int) = setindex!(parent(A), v, i)
Base.setindex!(A::AbstractConfig, v, I::Vararg{Int,N}) where N = setindex!(parent(A), v, I...)
Base.setindex!(A, X, I...) = setindex!(parent(A), X, I...)

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
