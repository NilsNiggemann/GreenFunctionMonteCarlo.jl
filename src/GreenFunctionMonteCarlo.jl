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

    include("AbstractTypes/AbstractTypes.jl")

    export AbstractWalkerEnsemble, AbstractPropagator, AbstractMove, AbstractConfig, AbstractHilbertSpace, AbstractOperator, AbstractGuidingFunction, AbstractGFMCProblem

    export ZeroDiagOperator

    export propagateWalkers!, fulfills_constraint, InverseMove, apply!

    include("Variational/EqualWeightSuperposition.jl")
    
    include("DefaultFormalism/BosonicConfig.jl")
    export BosonConfig, BosonHilbertSpace, OccupationNumberConstraint, HardCoreConstraint, fulfills_constraint

    include("DefaultFormalism/WalkerEnsemble.jl")
    
    include("DefaultFormalism/LocalOperator.jl")
    export LocalOperator, localOperator, SparseMove
    include("DefaultFormalism/ManyWalkerGFMC.jl")
    include("DefaultFormalism/ContinuousTimePropagator.jl")
    export ContinuousTimePropagator
end # module
