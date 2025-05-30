"""
    struct BasicAccumulator{T_high<:AbstractFloat} <: AbstractObserver

A structure that represents a basic accumulator for observing
energy and weights during simulations. Instead of storing observables at each step, the expectation values are computed on the fly, reducing the storage requirements significantly. Note that it is advisable to use a good guess of the average weight in the propagator to reduce numerical precision loss.

# Type Parameters
- `T_high`: The floating-point type used for observables, 
  constrained to subtypes of `AbstractFloat`.
"""
struct BasicAccumulator{T_high<:AbstractFloat} <: AbstractObserver
    TotalWeights::CircularArrays.CircularVector{T_high, Vector{T_high}}
    energies::CircularArrays.CircularVector{T_high, Vector{T_high}}
    Gnps::CircularArrays.CircularMatrix{T_high, Matrix{T_high}}
    reconfigurationTable::CircularArrays.CircularMatrix{Int, Matrix{Int}}
    PopulationMatrix::CircularArrays.CircularMatrix{Int, Matrix{Int}}
    en_numerator::Vector{T_high}
    Gnp_denominator::Vector{T_high}
end

function BasicAccumulator(filename,m_proj::Integer,NWalkers::Integer)
    p_proj = 2m_proj

    energies = CircularArrays.CircularArray(zeros(p_proj))
    TotalWeights = CircularArrays.CircularArray(zeros(p_proj))
    reconfigurationTable = CircularArrays.CircularArray(zeros(Int,NWalkers,p_proj))
    PopulationMatrix = CircularArrays.CircularArray(zeros(Int,NWalkers,m_proj))
    Gnps = CircularArrays.CircularArray(zeros(Float64,p_proj,p_proj))

    en_numerator = maybe_MMap_array(filename,"en_numerator",Float64,(m_proj,))
    Gnp_denominator = maybe_MMap_array(filename,"Gnp_denominator",Float64,(m_proj,))

    return BasicAccumulator(TotalWeights,energies,Gnps,reconfigurationTable,PopulationMatrix,en_numerator,Gnp_denominator)
end

BasicAccumulator(m_proj::Integer,NWalkers::Integer) = BasicAccumulator(nothing,m_proj,NWalkers)

function saveObservables_before!(Observables::BasicAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Hxx = get_diagonal(H)
    energies = Observables.energies
    TotalWeights = Observables.TotalWeights

    update_energies_TotalWeights!(energies,TotalWeights,i,Walkers,Hxx)

    Gnps = Observables.Gnps
    en_numerator = Observables.en_numerator
    Gnp_denominator = Observables.Gnp_denominator

    updateGnp!(Gnps,TotalWeights,i)
    NSites = length(getConfig(Walkers,1))
    getEnergy_step!(en_numerator,Gnp_denominator,Gnps,energies,i,NSites)
    return nothing
end

function getEnergy_step!(en_numerator::AbstractVector,Gnp_denominator::AbstractVector,Gnp::CircularArrays.CircularMatrix,localEnergies::AbstractVector,n::Integer,NSites::Integer)
    Nsites⁻¹ = 1/NSites
    for p in eachindex(en_numerator)
        n > p || continue
        en_numerator[p] += Gnp[n,p]*localEnergies[n]*Nsites⁻¹
        Gnp_denominator[p] += Gnp[n,p]
    end
    return en_numerator
end

function updateGnp!(Gnp,TotalWeights,n)
    nMax,pMax = size(Gnp)
    Gnp[n,1] = TotalWeights[n]>0 #zero projection order
    Gnp[n,2] = TotalWeights[n]

    for p in 3:pMax
        if n < p
            Gnp[n,p] = 0
            continue
        end
        Gnp[n,p] = Gnp[n-1,p-1]*Gnp[n,2]
    end
    return
end

function saveObservables_after!(Observables::BasicAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Observables.reconfigurationTable[:,i] .= get_reconfigurationList(reconfiguration)
    return nothing
end
