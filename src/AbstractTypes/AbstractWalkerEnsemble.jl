"""
    AbstractWalkerEnsemble

An abstract type that represents an ensemble of configurations (i.e. spins, bosons on a lattice for each walker) in the context of the Green Function Monte Carlo project. 

# Interface (required)
- getConfig(X::AbstractWalkerEnsemble,α):  get the configuration of the α-th walker in the ensemble.
- getMoveWeights(X::AbstractWalkerEnsemble,α): get the weights of the moves for the α-th walker in the ensemble.
- getBuffer(X::AbstractWalkerEnsemble,α): get the buffer associated with the α-th walker in the ensemble.
- getWalkerWeights(X::AbstractWalkerEnsemble): get the weights of the Walkers in the ensemble.
# Interface (optional)
"""
abstract type AbstractWalkerEnsemble end

function getConfig end
function getMoveWeights end
function getWalkerWeights end
function getBuffer end
# function getReconfigurationList end
NWalkers(X::AbstractWalkerEnsemble) = length(getConfigs(X))


