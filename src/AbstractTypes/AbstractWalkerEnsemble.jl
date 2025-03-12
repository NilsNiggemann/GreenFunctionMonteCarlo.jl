"""
    AbstractWalkerEnsemble

An abstract type that represents an ensemble of configurations (i.e. spins, bosons on a lattice for each walker) in the context of the Green Function Monte Carlo project. 

# Interface (required)
- eachConfig(X::AbstractWalkerEnsemble): iterate over the configurations of the ensemble.
- getWeights(X::AbstractWalkerEnsemble): get the weights of the configurations in the ensemble.
- getBuffers(X::AbstractWalkerEnsemble): get the buffers associated with the configurations in the ensemble.

# Interface (optional)
"""
abstract type AbstractWalkerEnsemble end

getConfigs(X::AbstractWalkerEnsemble) = throw(MethodError(getConfigs, (X,)))
getWeights(X::AbstractWalkerEnsemble) = throw(MethodError(getWeights, (X,)))
getBuffers(X::AbstractWalkerEnsemble) = throw(MethodError(getBuffers, (X,)))
NWalkers(X::AbstractWalkerEnsemble) = length(getConfigs(X))


