struct BosonHilbertSpace{Cons<:Tuple} <: AbstractHilbertSpace
    num_sites::Int
    Constraints::Cons
end

struct BosonConfig{T,N,Arr<:AbstractArray{T,N}} <: AbstractConfig{T,N}
    data::Arr
end

Base.parent(x::BosonConfig) = x.data

struct OccupationNumberConstraint <: AbstractConstraint
    min_occupation::Int
    max_occupation::Int
end

function (c::OccupationNumberConstraint)(x::BosonConfig)
    return all(c.min_occupation <= xi <= c.max_occupation for xi in x.data)
end

constraints(H::BosonHilbertSpace) = H.Constraints

function fulfills_constraint(x::BosonConfig, H::BosonHilbertSpace)
    return all(c(x) for c in constraints(H))
end