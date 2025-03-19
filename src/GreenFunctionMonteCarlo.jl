"""
Green function Monte Carlo is a method to sample from the ground state of sign-problem free Hamiltonians, i.e. those which can be writte as ```math H_{xx'} ≤ 0``` if ```math x != x'```.
"""
module GreenFunctionMonteCarlo

    import Random
    import SparseArrays
    import StaticArrays as SA
    import RecursiveArrayTools
    import StatsBase
    import ChunkSplitters
    import SmallCollections
    import HDF5
    import Statistics
    
    include("utils.jl")
    export createMMapArray, readMMapArray

    include("AbstractTypes/AbstractTypes.jl")
    
    export AbstractWalkerEnsemble, AbstractPropagator, AbstractMove, AbstractConfig, AbstractHilbertSpace, AbstractOperator, AbstractGuidingFunction, AbstractGFMCProblem, AbstractObserver, AbstractParallelizationScheme, AbstractConstraint, AbstractSignFreeOperator

    export ZeroDiagOperator

    export propagateWalkers!, fulfills_constraint, InverseMove, apply!

    export SingleThreaded, MultiThreaded
    
    include("Variational/EqualWeightSuperposition.jl")
    
    include("DefaultFormalism/BosonicConfig.jl")
    export BosonConfig, BosonHilbertSpace, OccupationNumberConstraint, HardCoreConstraint, fulfills_constraint

    include("DefaultFormalism/WalkerEnsemble.jl")
    
    include("DefaultFormalism/LocalOperator.jl")
    export LocalOperator, localOperator, SparseMove

    include("DefaultFormalism/MinimalReconfiguration.jl")
    include("DefaultFormalism/ConfigObserver.jl")
    export ConfigObserver
    include("DefaultFormalism/ManyWalkerGFMC.jl")
    export NoObserver, runGFMC!, GFMCProblem

    include("DefaultFormalism/ContinuousTimePropagator.jl")
    
    export ContinuousTimePropagator
end # module
