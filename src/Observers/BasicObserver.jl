
"""
    struct BasicObserver{DT<:AbstractFloat} <: AbstractObserver

Saves the energies, weights and reconfigurations of the walkers during the Monte Carlo simulation.

# Type Parameters
- `T`: Data type for the energies and weights.
"""
struct BasicObserver{T} <: AbstractObserver
    energies::Vector{T}
    TotalWeights::Vector{T}
    reconfigurationTable::Matrix{Int}
end
""" 
    function BasicObserver(filename,NSteps::Integer,NWalkers::Integer)
Creates a BasicObserver, storing energies, weights and reconfigurations. Creates a file under filename to store the data.
If no filename (or nothing) is provided, store the data in memory.
WARNING: Deleting a file that is being memory-mapped by a BasicObserver or any other observer will result in a crash when the observer is written to!
"""
function BasicObserver(filename,NSteps::Integer,NWalkers::Integer)
    energies = maybe_MMap_array(filename,"energies",Float64,(NSteps,))
    TotalWeights = maybe_MMap_array(filename,"TotalWeights",Float64,(NSteps,))
    reconfigurationTable = maybe_MMap_array(filename,"reconfigurationTable",Int,(NWalkers,NSteps))
    return BasicObserver(energies,TotalWeights,reconfigurationTable)
end
BasicObserver(NSteps::Integer,NWalkers::Integer) = BasicObserver(nothing,NSteps,NWalkers)

function saveObservables_before!(Observables::BasicObserver,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Hxx = get_diagonal(H)
    energies = Observables.energies
    TotalWeights = Observables.TotalWeights
    update_energies_TotalWeights!(energies,TotalWeights,i,Walkers,Hxx)
    return nothing
end

function update_energies_TotalWeights!(energies,TotalWeights,i,Walkers::AbstractWalkerEnsemble,Hxx::DiagonalOperator)
    WalkerWeights = getWalkerWeights(Walkers)
    energies[i] = getLocalEnergyWalkers_before(Walkers,Hxx)
    TotalWeights[i] = Statistics.mean(WalkerWeights)
    return nothing
end

function saveObservables_after!(Observables::BasicObserver,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Observables.reconfigurationTable[:,i] .= get_reconfigurationList(reconfiguration)
    return nothing
end
