"""
    struct WalkerAVGObserver{T} <: AbstractObserver

Saves the average of configurations over all walkers during the Monte Carlo simulation at each step.

# Type Parameters
- `T`: Data type for the configurations.
"""
struct WalkerAVGObserver{T} <: AbstractObserver
    average_configs::T
end

function WalkerAVGObserver(filename,config::AbstractConfig{T,N},NSteps::Integer;dtype=Float64) where {T,N}
    confSize = size(config)
    average_configs = maybe_MMap_array(filename,"average_configs",dtype,(confSize...,NSteps))
    return WalkerAVGObserver(average_configs)
end
WalkerAVGObserver(config::AbstractConfig,NSteps::Integer;kwargs...) = WalkerAVGObserver(nothing,config,NSteps;kwargs...)

function saveObservables_after!(Observables::WalkerAVGObserver,n,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    average_configs = Observables.average_configs
    for α in eachindex(Walkers)
        conf = parent(getConfig(Walkers,α))
        LoopVectorization.@turbo for i in eachindex(conf)
            average_configs[i,n] += conf[i]
        end
    end
    average_configs[:,n] ./= length(eachindex(Walkers))
    return nothing
end