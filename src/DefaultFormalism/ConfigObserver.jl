"""
    struct ConfigObserver{DT<:AbstractFloat,T} <: AbstractObserver

Saves the configurations, energies, and weights of the walkers during the Monte Carlo simulation.

# Type Parameters
- `DT<:AbstractFloat`: The floating-point type used for numerical computations (e.g., `Float64` or `Float32`).
- `T`: A generic type parameter representing additional data or configuration associated with the observer.

# Usage
This type is typically used internally within the framework to monitor and process configurations during the Monte Carlo simulation.
"""
struct ConfigObserver{DT<:AbstractFloat,T} <: AbstractObserver
    energies::Vector{DT}
    SaveConfigs::T
    TotalWeights::Vector{DT}
    reconfigurationTable::Matrix{Int}
end


function ConfigObserver(filename,config::AbstractConfig{T,N},NSteps::Integer,NWalkers::Integer) where {T,N}
    confSize = size(config)
    energies = maybe_MMap_array(filename,"energies",Float64,(NSteps,))
    TotalWeights = maybe_MMap_array(filename,"TotalWeights",Float64,(NSteps,))
    reconfigurationTable = maybe_MMap_array(filename,"reconfigurationTable",Int,(NWalkers,NSteps))
    SaveConfigs = maybe_MMap_array(filename,"SaveConfigs",T,(confSize...,NWalkers,NSteps))
    return ConfigObserver(energies,SaveConfigs,TotalWeights,reconfigurationTable)
end
ConfigObserver(config::AbstractConfig,NSteps::Integer,NWalkers::Integer) = ConfigObserver(nothing,config,NSteps,NWalkers)

get_reconfigurationTable(Observables::ConfigObserver) = Observables.reconfigurationTable
get_energies(Observables::ConfigObserver) = Observables.energies
get_TotalWeights(Observables::ConfigObserver) = Observables.TotalWeights

function update_energies_TotalWeights!(energies,TotalWeights,i,Walkers::AbstractWalkerEnsemble,Hxx::DiagonalOperator)
    WalkerWeights = getWalkerWeights(Walkers)
    energies[i] = getLocalEnergyWalkers_before(Walkers,Hxx)
    TotalWeights[i] = Statistics.mean(WalkerWeights)
    return nothing
end

function saveObservables_before!(Observables::ConfigObserver,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Hxx = get_diagonal(H)
    energies = get_energies(Observables)
    TotalWeights = get_TotalWeights(Observables)
    update_energies_TotalWeights!(energies,TotalWeights,i,Walkers,Hxx)
    return nothing
end

function saveObservables_after!(Observables::ConfigObserver,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    SaveConfigs = Observables.SaveConfigs

    for α in eachindex(Walkers)
        SaveConfigs[:,α,i] .= getConfig(Walkers,α)
    end
    Observables.reconfigurationTable[:,i] .= get_reconfigurationList(reconfiguration)
    return nothing
end
