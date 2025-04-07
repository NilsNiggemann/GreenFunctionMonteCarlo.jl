mutable struct EnsembleLogger_2{T} <: AbstractLogger
    obs_Loggers::Vector{SilentLogger_1}
    p::Union{Nothing,ProgressMeter.Progress}
    options::T
end
function getIter(logger::EnsembleLogger_2)
    return sum(log.current_step for log in logger.obs_Loggers)
end
function EnsembleLogger_2(;kwargs...)
    loggers = SilentLogger_1[]
    return EnsembleLogger_2(loggers,nothing,kwargs)
end

function setup_logger!(logger::EnsembleLogger_2, Walkers::AbstractWalkerEnsemble, Observables::Any, reconfiguration::AbstractReconfigurationScheme)
    logger.p = ProgressMeter.Progress(length(range)*length(P.problems),dt=0.01;output = stderr,showspeed=true, enabled = !is_logging(stderr),logger.options...)
    logger.obs_Loggers = [SilentLogger_1(0) for _ in P.problems]
    return logger
end

function write_log(logger::EnsembleLogger_2,i, range, Walkers::AbstractWalkerEnsemble, Observables::Any, reconfiguration::AbstractReconfigurationScheme)
    # showValFunc() = generate_showvalues(i, Walkers, Observables, reconfiguration)
    # iteration = getIter(logger)
    # println("Iteration: $iteration")
    # ProgressMeter.next!(logger.p;desc = "GFMC Ensemble step... $iteration/$(length(range)*length(P.problems))")
    ProgressMeter.update!(logger.p,i;desc = "GFMC Ensemble step... $iteration/$(length(range)*length(P.problems))")
    try
        1
    catch e
        @warn e maxlog=1
    end
    
end