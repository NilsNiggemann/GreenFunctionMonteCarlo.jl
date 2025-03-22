"""
    NoLogger <: AbstractLogger

A placeholder implementation that does not perform any logging.
# Example
"""
struct NoLogger <: AbstractLogger end
write_log(logger::NoLogger,i,range,Walkers,Observables,reconfiguration) = nothing

"""
    struct SimpleLogger <: AbstractLogger

A simple logger implementation that inherits from `AbstractLogger`.

# Fields
- `n_report::Int`: The number of steps between each log report.

This logger can be used to control and manage logging behavior in a straightforward manner.
"""
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

function logstring_energies(Observables::AbstractObserver,i)
    return ""
end

function write_log(logger::SimpleLogger,i,range,Walkers::AbstractWalkerEnsemble,Observables::Any,reconfiguration::AbstractReconfigurationScheme)
    if i % logger.n_report == 0
        try
            showvals = generate_showvalues(i,Walkers,Observables,reconfiguration)
            for (name, value) in showvals
                (name,value) === _empty_log() && continue
                print("$name: $value  ")
            end
            println()
        catch e
            @warn e maxlog=1
        end
    end
end
