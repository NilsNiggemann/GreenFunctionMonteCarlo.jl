"""
    struct ConfigurationObserver{T} <: AbstractObserver

Saves the configurations of the walkers during the Monte Carlo simulation.

# Type Parameters
- `T`: Data type for the configurations.
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

function ConfigObserver(filename,config::AbstractConfig{T,N},NSteps::Integer,NWalkers::Integer) where {T,N}
    return CombinedObserver((;BasicObserver = BasicObserver(filename,NSteps,NWalkers),ConfigurationObserver = ConfigurationObserver(filename,config,NSteps,NWalkers)))
end
ConfigObserver(config::AbstractConfig,NSteps::Integer,NWalkers::Integer) = ConfigObserver(nothing,config,NSteps,NWalkers)