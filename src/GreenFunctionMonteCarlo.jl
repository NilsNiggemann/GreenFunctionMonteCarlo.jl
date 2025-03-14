"""
Green function Monte Carlo is a method to sample from the ground state of sign-problem free Hamiltonians, i.e. those which can be writte as ```math H_{xx'} ≤ 0``` if ```math x != x'```.
"""
module GreenFunctionMonteCarlo

    import SparseArrays
    import StaticArrays as SA
    import RecursiveArrayTools
    import StatsBase

    include("AbstractTypes/AbstractTypes.jl")

    export AbstractWalkerEnsemble, AbstractPropagator, AbstractMove, AbstractConfig, AbstractHilbertSpace, AbstractOperator, AbstractGuidingFunction, AbstractGFMCProblem
    export propagateWalkers!, fulfills_constraints, InverseMove, apply!

    include("Variational/EqualWeightSuperposition.jl")
    
    include("DefaultFormalism/BosonicConfig.jl")
    export BosonConfig, BosonHilbertSpace, OccupationNumberConstraint, fulfills_constraint

    include("DefaultFormalism/LocalOperator.jl")
    include("DefaultFormalism/ManyWalkerGFMC.jl")
    include("DefaultFormalism/ContinuousTimePropagator.jl")
end # module
