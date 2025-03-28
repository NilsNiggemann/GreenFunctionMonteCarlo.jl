"""
    struct BosonHilbertSpace{Cons<:AbstractConstraint} <: AbstractHilbertSpace

A structure representing the Hilbert space for bosonic systems.

# Type Parameters
- `Cons<:AbstractConstraint`: A type that defines the constraints applied to the bosonic Hilbert space. This allows for flexible customization of the space's properties.

# Supertype
- `AbstractHilbertSpace`: This structure is a subtype of `AbstractHilbertSpace`, indicating that it represents a specific type of quantum mechanical Hilbert space.
"""
struct BosonHilbertSpace{Cons<:AbstractConstraint} <: AbstractHilbertSpace
    num_sites::Int
    Constraint::Cons
end
BosonHilbertSpace(num_sites::Int) = BosonHilbertSpace(num_sites, NoConstraint())

constraint(H::BosonHilbertSpace) = H.Constraint

"""
    struct BosonConfig{T, N, Arr<:AbstractArray{T, N}} <: AbstractConfig{T, N}

Represents a configuration of bosonic particles in a given formalism. This structure is parameterized by:

- `T`: The type of the elements in the configuration (e.g., `Bool` or `UInt8`).
- `N`: The dimensionality of the configuration space.
- `Arr<:AbstractArray{T, N}`: The type of the array used to store the configuration data, which must be a subtype of `AbstractArray` with element type `T` and dimensionality `N`.

This type is a subtype of `AbstractConfig{T, N}`, which provides a common interface for configurations in the system and implements the AbstractArray interface.
"""
struct BosonConfig{T,N,Arr<:AbstractArray{T,N}} <: AbstractConfig{T,N}
    data::Arr
end
Base.parent(x::BosonConfig) = x.data
Base.copy(x::BosonConfig) = BosonConfig(copy(x.data))
Base.size(H::BosonHilbertSpace) = (H.num_sites,)

"""
    HardCoreConstraint <: AbstractConstraint

A struct representing a hard-core constraint in the system. This constraint 
enforces that certain configurations are not allowed, typically used in 
bosonic systems to model hard-core interactions where particles cannot 
occupy the same state or position.

This type is a subtype of `AbstractConstraint`, which serves as a base 
for defining various constraints in the system.
Useful for simulating spin-1/2 systems.
"""
struct HardCoreConstraint <: AbstractConstraint end
struct OccupationNumberConstraint <: AbstractConstraint
    min_occupation::Int
    max_occupation::Int
end

get_boson_datatype(C::AbstractConstraint) = UInt8
get_boson_datatype(C::HardCoreConstraint) = Bool

function BosonConfig(H::BosonHilbertSpace)
    C = constraint(H)
    BosonConfig(zeros(get_boson_datatype(C), size(H)))
end

function (c::OccupationNumberConstraint)(x::BosonConfig)
    return all(c.min_occupation <= xi <= c.max_occupation for xi in x)
end
(c::HardCoreConstraint)(x::BosonConfig{Bool}) = true
(c::HardCoreConstraint)(x::BosonConfig) = all(0 <= xi <= 1 for xi in x)