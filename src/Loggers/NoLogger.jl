"""
    NoLogger <: AbstractLogger

A Logger that suppresses all logging output. 
"""
struct NoLogger <: AbstractLogger end
write_log(logger::NoLogger,args...;kwargs...) = nothing
write_log(logger::Nothing,args...;kwargs...) = nothing
