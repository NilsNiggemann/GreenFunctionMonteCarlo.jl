empty_log() = ("empty", "empty")
is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
default_logger() = is_logging(stderr) ? SimpleLogger(10) : ProgressBarLogger()

log_observable(O::AbstractObserver,i) = (empty_log(),)

function log_walker_survival_ratio(reconfigurationTable,i)
    RL = @view reconfigurationTable[:,i]
    survivingWalkers = unique(RL)
    ratio = length(survivingWalkers)/length(RL)
    return ("Walker survival rate",strd(ratio))
end

log_observable(O::BasicObserver,i) = (log_walker_survival_ratio(O.reconfigurationTable,i),log_Obs_energy(O,i),log_Obs_weights(O,i))
log_Obs_energy(O::BasicObserver,i) = ("eloc",strd(O.energies[i]))
log_Obs_weights(O::BasicObserver,i) = ("w_avg",strd(O.TotalWeights[i]))

function log_observable(O::CombinedObserver,i)
    vals = [empty_log()]
    for obs in O.Observers
        appendlog!(vals,log_observable(obs,i))
    end
    return vals
end

appendlog!(log::Vector{Tuple{String,String}},(a,b)::Tuple{A,B}) where {A,B} = push!(log,(string(a),string(b)))

function appendlog!(log::Vector{Tuple{String,String}},vals)
    for val in vals
        appendlog!(log,val)
    end
    return vals
end

function generate_showvalues(i, Walkers::AbstractWalkerEnsemble,Observer::AbstractObserver,reconfiguration::AbstractReconfigurationScheme)
    vals = log_observable(Observer,i)
    vals = filterObs!(vals)
    return vals
end
function filterObs!(vals::Any)
    return filter(x -> x != empty_log(),vals)
end
# function filterObs!(vals::Vector)
#     return filter!(x -> x != empty_log(),vals)
# end
function filterObs!(vals::Tuple{String,String})
    vals == empty_log() && return ()
    return vals
end