struct NoLogger <: AbstractLogger end
write_log(logger::NoLogger,i,range,Walkers,Observables,reconfiguration) = nothing

struct SimpleLogger <: AbstractLogger
    n_report::Int
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
