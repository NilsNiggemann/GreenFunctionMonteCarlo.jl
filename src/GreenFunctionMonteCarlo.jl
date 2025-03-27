module GreenFunctionMonteCarlo
    # Use README as the docstring of the module:
    @doc read(joinpath(dirname(@__DIR__), "README.md"), String) GreenFunctionMonteCarlo
    
    import Random
    import SparseArrays
    import StaticArrays as SA
    import RecursiveArrayTools
    import StatsBase
    import ChunkSplitters
    import SmallCollections
    import HDF5
    import Statistics
    import LinearAlgebra
    import LoopVectorization
    import ProgressMeter
    import FFTW

    include("utils.jl")
    export createMMapArray, readMMapArray

    include("AbstractTypes/AbstractTypes.jl")

    export AbstractWalkerEnsemble, AbstractPropagator, AbstractMove, AbstractConfig, AbstractHilbertSpace, AbstractOperator, AbstractGuidingFunction, AbstractGFMCProblem, AbstractObserver, AbstractParallelizationScheme, AbstractConstraint, AbstractSignFreeOperator, AbstractLogger

    export ZeroDiagOperator

    export propagateWalkers!, fulfills_constraint, InverseMove, apply!, get_params

    export SingleThreaded, MultiThreaded

    include("Observers/BasicObserver.jl")
    include("Observers/CombinedObserver.jl")
    include("Observers/ConfigObserver.jl")
    export BasicObserver, ConfigurationObserver, ConfigObserver

    include("DefaultFormalism/BosonicConfig.jl")
    export BosonConfig, BosonHilbertSpace, OccupationNumberConstraint, HardCoreConstraint, fulfills_constraint

    include("DefaultFormalism/WalkerEnsemble.jl")

    include("DefaultFormalism/LocalOperator.jl")
    export LocalOperator, DiagOperator, localOperator, SparseMove

    include("Reconfiguration/MinimalReconfiguration.jl")

    include("DefaultFormalism/ManyWalkerGFMC.jl")
    export NoObserver, runGFMC!, GFMCProblem, ProblemEnsemble

    include("DefaultFormalism/ContinuousTimePropagator.jl")
    export ContinuousTimePropagator

    include("DefaultFormalism/Evaluation.jl")
    export getEnergies, precomputeNormalizedAccWeight

    include("Variational/EqualWeightSuperposition.jl")
    export EqualWeightSuperposition

    include("Variational/Jastrow.jl")
    export Jastrow

    include("Variational/NaiveFunction.jl")
    export NaiveFunction


    include("Loggers/LoggerUtils.jl")

    include("Loggers/NoLogger.jl")
    export NoLogger

    include("Loggers/SimpleLogger.jl")
    export SimpleLogger

    include("Loggers/ProgressBarLogger.jl")
    export ProgressBarLogger

    include("Observables/computeObservables.jl")
    include("Observables/CorrelationFunction.jl")
end # module
