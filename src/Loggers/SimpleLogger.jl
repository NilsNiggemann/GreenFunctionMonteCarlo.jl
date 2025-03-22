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
