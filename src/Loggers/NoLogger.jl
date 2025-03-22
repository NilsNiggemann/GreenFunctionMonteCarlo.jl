"""
    NoLogger <: AbstractLogger

A placeholder implementation that does not perform any logging.
# Example
"""
struct NoLogger <: AbstractLogger end
write_log(logger::NoLogger,i,range,Walkers,Observables,reconfiguration) = nothing
