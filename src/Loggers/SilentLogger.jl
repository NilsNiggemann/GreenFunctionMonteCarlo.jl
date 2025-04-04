mutable struct SilentLogger_1 <: AbstractLogger
    current_step::Int
    showvalues::Vector{Tuple{String, String}}
end

function write_log(logger::SilentLogger_1, i, range, Walkers::AbstractWalkerEnsemble, Observables::Any, reconfiguration::AbstractReconfigurationScheme)
    showvals = generate_showvalues(i, Walkers, Observables, reconfiguration)
    logger.current_step = i
    logger.showvalues = showvals
end
