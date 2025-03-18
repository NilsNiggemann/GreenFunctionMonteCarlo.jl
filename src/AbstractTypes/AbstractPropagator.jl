"""
    AbstractPropagator

An abstract type representing a propagator. Examples include discrete or continuous time methods.

# Interface:
- `propagateWalkers!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, propagator::AbstractPropagator, parallelization::AbstractParallelizationScheme, RNG::Random.AbstractRNG = Random.default_rng())`: Propagate the walkers in the ensemble using the specified moves and parameters.
"""
abstract type AbstractPropagator end

# propagateWalkers!(X::AbstractWalkerEnsemble, moves, P::AbstractPropagator,params) = throw(MethodError(propagateWalkers!, (X, moves, P, params)))

""" Propagates walker ensemble `X` using a collection of moves according to the rules specified by the propagator `P`. 

# Usage:    
`propagateWalkers!(WE::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, propagator::AbstractPropagator, parallelization::AbstractParallelizationScheme, RNG::Random.AbstractRNG = Random.default_rng())`
"""
function propagateWalkers! end