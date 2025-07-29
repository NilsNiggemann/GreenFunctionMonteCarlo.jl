"""
    struct WalkerAVGObserver{T} <: AbstractObserver

Saves the average of configurations over all walkers during the Monte Carlo simulation at each step.

# Type Parameters
- `T`: Data type for the configurations.
"""
struct WalkerAVGObserver{T} <: AbstractObserver
    average_configs::T
    bin_elements::Int
end

function WalkerAVGObserver(filename,config::AbstractConfig{T,N},NSteps::Integer;dtype=Float64,bin_elements=1) where {T,N}
    confSize = size(config)
    average_configs = maybe_MMap_array(filename,"average_configs",dtype,(confSize...,NSteps÷bin_elements))
    return WalkerAVGObserver(average_configs, bin_elements)
end
WalkerAVGObserver(config::AbstractConfig,NSteps::Integer;kwargs...) = WalkerAVGObserver(nothing,config,NSteps;kwargs...)

function saveObservables_after!(Observables::WalkerAVGObserver,n,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    average_configs = Observables.average_configs
    Nbins = last(size(average_configs))
    n_idx = get_bin_index(n,Nbins,Observables.bin_elements)
    for α in eachindex(Walkers)
        conf = parent(getConfig(Walkers,α))
        LoopVectorization.@turbo for i in eachindex(conf)
            average_configs[i,n_idx] += conf[i]
        end
    end
    average_configs[:,n_idx] ./= length(eachindex(Walkers))
    return nothing
end