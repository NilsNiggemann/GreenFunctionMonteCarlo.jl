"""
    struct ConfigurationObserver{T} <: AbstractObserver

Saves the configurations of the walkers during the Monte Carlo simulation.

# Type Parameters
- `T`: Data type for the configurations.

# See Also
- [`ObservableAccumulator`](@ref)
"""
struct ConfigurationObserver{T} <: AbstractObserver
    SaveConfigs::T
end

function ConfigurationObserver(filename,config::AbstractConfig{T,N},NSteps::Integer,NWalkers::Integer) where {T,N}
    confSize = size(config)
    SaveConfigs = maybe_MMap_array(filename,"SaveConfigs",T,(confSize...,NWalkers,NSteps))
    return ConfigurationObserver(SaveConfigs)
end
ConfigurationObserver(config::AbstractConfig,NSteps::Integer,NWalkers::Integer) = ConfigurationObserver(nothing,config,NSteps,NWalkers)

function saveObservables_after!(Observables::ConfigurationObserver,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    SaveConfigs = Observables.SaveConfigs
    for α in eachindex(Walkers)
        SaveConfigs[:,α,i] .= getConfig(Walkers,α)
    end
    return nothing
end

"""
    ConfigObserver(filename, config::AbstractConfig{T,N}, NSteps::Integer, NWalkers::Integer) where {T,N}

Creates a combined observer that tracks energy, weight, reconfigurations and configurations of the walkers during the Monte Carlo simulation. The observer saves the data to a file with the given filename.

# Arguments
- `filename::String`: The name of the file where the configuration data will be saved.
- `config::AbstractConfig{T,N}`: The initial configuration object representing the system state.
- `NSteps::Integer`: The number of simulation steps to be observed.
- `NWalkers::Integer`: The number of walkers (or particles) in the simulation.
"""
function ConfigObserver(filename,config::AbstractConfig{T,N},NSteps::Integer,NWalkers::Integer) where {T,N}
    return CombinedObserver((;BasicObserver = BasicObserver(filename,NSteps,NWalkers),ConfigurationObserver = ConfigurationObserver(filename,config,NSteps,NWalkers)))
end
ConfigObserver(config::AbstractConfig,NSteps::Integer,NWalkers::Integer) = ConfigObserver(nothing,config,NSteps,NWalkers)