using ProgressMeter
import GreenFunctionMonteCarlo as GFMC

struct ProgressBarLogger{T} <: AbstractLogger
    p::Ref{Union{Nothing,Progress}}
    dt::Float64
    options::T
end

"""
    ProgressBarLogger(; dt=0.1, kwargs...)

Creates a progress bar logger with customizable update intervals and additional options.

# Arguments
- `dt::Float64`: The time interval (in seconds) between progress bar updates. Defaults to `0.1`.
- `kwargs...`: Additional keyword arguments to customize the behavior of the progress bar logger. See `GreenFunctionMonteCarlo.ProgressMeter.Progress` for more details.

# Returns
A progress bar logger instance that can be used to track and display progress in a task.

# Example
"""
function ProgressBarLogger(;dt=0.1,kwargs...) 
    return ProgressBarLogger{typeof(kwargs)}(nothing, dt,kwargs)
end

is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")

function write_log(logger::ProgressBarLogger, i, range, Walkers::AbstractWalkerEnsemble, Observables::Any, reconfiguration)

    if isnothing(logger.p[]) # if logger.p is nothing, create a new progress bar
        logger.p[] = Progress(length(range),dt=logger.dt;output = stderr,showspeed=true, enabled = !is_logging(stderr),desc="running GFMC...",logger.options...)
    end
    
    next!(logger.p[],showvalues = () -> generate_showvalues(i, Walkers,Observables,reconfiguration))
end

function write_log(logger::ProgressBarLogger, i, range, Walkers::AbstractWalkerEnsemble, Observables::NoObserver, reconfiguration)
    if isnothing(logger.p[])
        logger.p[] = Progress(length(range),dt=logger.dt)
    end
    RL = get_reconfigurationList(reconfiguration)
    SurvWalkers = length(unique(RL))
    NW = NWalkers(Walkers)
    next!(logger.p[],showvalues = [("Walker Survival",SurvWalkers/NW)])
end