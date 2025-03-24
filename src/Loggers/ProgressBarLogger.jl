mutable struct ProgressBarLogger{T} <: AbstractLogger
    p::Union{Nothing,ProgressMeter.Progress}
    options::T
end

"""
    ProgressBarLogger(; kwargs...)

Creates a progress bar logger with customizable update intervals and additional options.

# Arguments
- `kwargs...`: Additional keyword arguments to customize the behavior of the progress bar logger. See `GreenFunctionMonteCarlo.ProgressMeter.Progress` for more details.

# Returns
A progress bar logger instance that can be used to track and display progress in a task.

# Example
"""
function ProgressBarLogger(;kwargs...)
    return ProgressBarLogger(nothing,kwargs)
end

function write_log(logger::ProgressBarLogger, i, range, Walkers::AbstractWalkerEnsemble, Observables::Any, reconfiguration)
    if isnothing(logger.p) || i == first(range) # initialize progress bar
        logger.p = ProgressMeter.Progress(length(range),dt=0.1;output = stderr,showspeed=true, enabled = !is_logging(stderr),logger.options...)
    end

    showValFunc() = generate_showvalues(i, Walkers, Observables, reconfiguration)

    ProgressMeter.next!(logger.p,showvalues = showValFunc();desc = "running GFMC... $i/$(last(range))")
    try
        1
    catch e
        @warn e maxlog=1
    end
    
end