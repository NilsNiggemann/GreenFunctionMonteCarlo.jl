# module ProgressBarLoggerExt
    # using GreenFunctionMonteCarlo, ProgressBars
    using ProgressMeter
    import GreenFunctionMonteCarlo as GFMC

    struct ProgressBarLogger <: AbstractLogger
        p::Ref{Union{Nothing,Progress}}
        dt::Float64
    end
    ProgressBarLogger(;dt=0.1) = ProgressBarLogger(nothing, dt)

    function write_log(logger::ProgressBarLogger, i, range, Walkers::AbstractWalkerEnsemble, Observables::Any, reconfiguration)
        if isnothing(logger.p[])
            logger.p[] = Progress(length(range),dt=logger.dt)
        end
        next!(logger.p[],showvalues = [("eloc",Observables.energies[i]), ("w_avg",Observables.TotalWeights[i])])
    end
# end
