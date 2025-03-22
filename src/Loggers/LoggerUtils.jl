
# _empty_log() = nothing
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
    vals = (("Iteration",i),log_walker_survival_ratio(reconfiguration,i),log_Obs_energy(Observer,i),log_Obs_weights(Observer,i))

    # return Iterators.filter(x -> x != _empty_log(),vals)
    # return Iterators.filter(x -> x != _empty_log(),vals)
end