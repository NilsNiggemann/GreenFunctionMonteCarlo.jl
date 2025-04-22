@doc raw"""
## Overview

`GreenFunctionMonteCarlo.jl` is a Julia package designed for performing Green Function Monte Carlo (GFMC) simulations on lattice models.

Presently, this package treats only Hamiltonians that are free of the sign problem, i.e. whose elements satisfy
```math

H_{x, x'} \leq 0 \quad \forall x \neq x'
```
where $H_{x, x'}$ is the matrix element of the Hamiltonian between two configurations (spins or bosons) $x$ and $x'$. 

## Usage: 
```
using GreenFunctionMonteCarlo, LinearAlgebra
NSites = 3
Nwalkers = 10
Hilbert = BosonHilbertSpace(NSites, HardCoreConstraint())
moves = eachcol(Bool.(I(NSites))) # each move flips a single spin
offdiagElements = -ones(NSites)
H = localOperator(moves, offdiagElements, DiagOperator(x->0), Hilbert)

prob = GFMCProblem(BosonConfig(Hilbert), Nwalkers, ContinuousTimePropagator(0.1); logψ = EqualWeightSuperposition(), H, Hilbert)
Observer = ConfigObserver(startConfig, NSteps, NWalkers) # Observer to measure the energy and configurations
runGFMC!(problem, NoObserver(), NStepsEquil) #run for NStepsEquil steps without observing to equilibrate
runGFMC!(problem, Observer, NSteps) #run for NSteps steps
```
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
    import LinearAlgebra
    import LoopVectorization
    import ProgressMeter
    import CircularArrays
    include("utils.jl")
    export createMMapArray, readMMapArray

    include("AbstractTypes/AbstractTypes.jl")

    export AbstractWalkerEnsemble, AbstractPropagator, AbstractMove, AbstractConfig, AbstractHilbertSpace, AbstractOperator, AbstractGuidingFunction, AbstractGFMCProblem, AbstractObserver, AbstractObservable, AbstractParallelizationScheme, AbstractConstraint, AbstractSignFreeOperator, AbstractLogger

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
    export LocalOperator, ZeroDiagOperator, DiagOperator, OneBodyDiagOperator, TwoBodyDiagOperator, DiagOperatorSum, localOperator, SparseMove

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
    export getObs_diagonal

    include("Observers/BasicAccumulator.jl")
    export BasicAccumulator
end # module
