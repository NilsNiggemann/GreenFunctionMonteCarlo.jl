struct NoLogger <: AbstractLogger end
write_log(logger::NoLogger,i,range,Walkers,Observables,reconfiguration) = nothing

struct SimpleLogger <: AbstractLogger
    n_report::Int
end

_empty_log() = ("", nothing)

function log_walker_survival_ratio(reconfiguration::MinimalReconfiguration,i)
    RL = get_reconfigurationList(reconfiguration)
    survivingWalkers = unique(RL)
    ratio = length(survivingWalkers)/length(RL)
    return ("Walker survival rate",ratio)
end
log_walker_survival_ratio(reconfiguration::AbstractReconfigurationScheme,i) = _empty_log()

log_Obs_energy(O::ConfigObserver,i) = ("eloc",O.energies[i])
log_Obs_energy(Observables::AbstractObserver,i) = _empty_log()

log_Obs_weights(O::ConfigObserver,i) = ("w_avg",O.TotalWeights[i])
log_Obs_weights(Observables::AbstractObserver,i) = _empty_log()

function generate_showvalues(i, Walkers::AbstractWalkerEnsemble,Observer::AbstractObserver,reconfiguration::AbstractReconfigurationScheme)
    return (("Iteration",i),log_walker_survival_ratio(reconfiguration,i),log_Obs_energy(Observer,i),log_Obs_weights(Observer,i))
end

function logstring_reconfiguration(Walkers::AbstractWalkerEnsemble,reconfiguration::AbstractReconfigurationScheme)
    RL = get_reconfigurationList(reconfiguration)
    survivingWalkers = unique(RL)
    str = "surviving Walkers: $(length(survivingWalkers)) out of $(NWalkers(Walkers))"
    return str
end

function logstring_energies(Observables::ConfigObserver,i)
    (;energies,TotalWeights) = Observables
    # en_mean = Statistics.mean(Observables.energies)
    # en_std = Statistics.std(Observables.energies)
    # TotalWeights_mean = Statistics.mean(Observables.TotalWeights)
    # str = "Energy: $(strd(en_mean)) ± $(strd(en_std)) w_avg: $(strd(TotalWeights_mean))"
    str = "e_local: $(strd(energies[i])) w_avg: $(strd(TotalWeights[i]))"
    return str
end

function logstring_energies(Observables::AbstractObserver,i)
    return ""
end

function write_log(logger::SimpleLogger,i,range,Walkers::AbstractWalkerEnsemble,Observables::Any,reconfiguration::AbstractReconfigurationScheme)
    if i % logger.n_report == 0
        try
            print("Iteration: $i of $(last(range))\t")
            print(logstring_reconfiguration(Walkers,reconfiguration),"\t")
            print(logstring_energies(Observables,i),"\n")
        catch e
            @warn e
        end
    end
end
