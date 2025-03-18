
"""
    AbstractReconfigurationScheme

An abstract type that serves as a base for defining different reconfiguration strategies 
in the context of the Green Function Monte Carlo framework. 
# Interface 
- `reconfigurateWalkers!(Walkers::AbstractWalkerEnsemble,reconfiguration::AbstractReconfigurationScheme,rng::Random.AbstractRNG)`
"""
abstract type AbstractReconfigurationScheme end