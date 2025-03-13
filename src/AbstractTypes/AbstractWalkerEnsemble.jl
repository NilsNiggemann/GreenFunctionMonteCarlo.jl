"""
    AbstractWalkerEnsemble

An abstract type that represents an ensemble of configurations (i.e. spins, bosons on a lattice for each walker) in the context of the Green Function Monte Carlo project. 

# Interface (required)
- eachConfig(X::AbstractWalkerEnsemble): iterate over the configurations of the ensemble.
- getWalkerWeights(X::AbstractWalkerEnsemble): get the weights of the configurations in the ensemble.
- getBuffers(X::AbstractWalkerEnsemble): get the buffers associated with the configurations in the ensemble.
- getReconfigurationList(X::AbstractWalkerEnsemble): get the reconfiguration list associated with the ensemble.
# Interface (optional)
"""
abstract type AbstractWalkerEnsemble end

function getConfigs end
function getWalkerWeights end
function getWeightLists end
function getBuffers end
function getReconfigurationList end
NWalkers(X::AbstractWalkerEnsemble) = length(getConfigs(X))


