struct BosonHilbertSpace{Cons<:AbstractConstraint} <: AbstractHilbertSpace
    num_sites::Int
    Constraint::Cons
end
BosonHilbertSpace(num_sites::Int) = BosonHilbertSpace(num_sites, NoConstraint())

constraint(H::BosonHilbertSpace) = H.Constraint
struct BosonConfig{T,N,Arr<:AbstractArray{T,N}} <: AbstractConfig{T,N}
    data::Arr
end
Base.parent(x::BosonConfig) = x.data
Base.copy(x::BosonConfig) = BosonConfig(copy(x.data))
Base.size(H::BosonHilbertSpace) = (H.num_sites,)

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