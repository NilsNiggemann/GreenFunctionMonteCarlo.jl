struct ConfigSaver{DT<:AbstractFloat,T,T2} <: AbstractObservables
    energies::Vector{DT}
    SaveConfigs::T
    TotalWeights::Vector{DT}
    reconfigurationTable::Matrix{Int}
end
get_reconfigurationTable(Observables::ConfigSaver) = Observables.reconfigurationTable
get_energies(Observables::ConfigSaver) = Observables.energies
get_TotalWeights(Observables::ConfigSaver) = Observables.TotalWeights

function update_energies_TotalWeights!(energies,TotalWeights,i,Walkers::AbstractWalkerEnsemble,propagator::AbstractPropagator)
    WalkerWeights = getWalkerWeights(Walkers)
    energies[i] = getLocalEnergyWalkers_before(Walkers,propagator::AbstractPropagator)
    TotalWeights[i] = mean(WalkerWeights)
    return nothing
end

function saveObservables_before!(Observables::ConfigSaver,i,Walkers::AbstractWalkerEnsemble,propagator::AbstractPropagator)
    energies = get_energies(Observables)
    TotalWeights = get_TotalWeights(Observables)
    update_energies_TotalWeights!(energies,TotalWeights,i,Walkers,propagator::AbstractPropagator)
    return nothing
end

function saveObservables_after!(Observables::ConfigSaver,i,Walkers::AbstractWalkerEnsemble)
    SaveConfigs = Observables.SaveConfigs

    for α in eachindex(Walkers)
        SaveConfigs[:,α,i] .= getConfig(Config,α)
    end
    return nothing
end