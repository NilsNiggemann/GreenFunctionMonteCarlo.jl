struct EnsembleLogger_1 <: AbstractLogger
    obs_Loggers::Vector{SilentLogger_1}
    p::ProgressMeter.Progress
end
EnsembleLogger_1()
function write_log(logger::EnsembleLogger_1, i, range, Walkers::AbstractWalkerEnsemble, Observables::Any, reconfiguration)
    if isnothing(logger.p) || i == first(range) # initialize progress bar
        logger.p = ProgressMeter.Progress(length(range),dt=0.1;output = stderr,showspeed=true, enabled = !is_logging(stderr),logger.options...)
    end

    showValFunc() = generate_showvalues(i, Walkers, Observables, reconfiguration)

    ProgressMeter.next!(logger.p,showvalues = showValFunc();desc = "GFMC step... $i/$(last(range))")
    try
        1
    catch e
        @warn e maxlog=1
    end
    
end