"""
    CombinedObserver{T}

A struct that represents a combined observer, which is a collection of multiple observers grouped together. 
This allows for observing multiple quantities or behaviors simultaneously.

# Fields
- `Observers::T`: A collection containing the individual observers.

# Usage
This struct is useful for combining multiple observer objects into a single entity, enabling them to be managed and used collectively.
"""
struct CombinedObserver{T} <: AbstractObserver
    Observers::T
end

function saveObservables_before!(Observer::CombinedObserver,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    for obs in Observer.Observers
        saveObservables_before!(obs,i,Walkers,H,reconfiguration)
    end
    return nothing
end

function saveObservables_after!(Observer::CombinedObserver,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    for obs in Observer.Observers
        saveObservables_after!(obs,i,Walkers,H,reconfiguration)
    end
    return nothing
end

function getObs(Observer::CombinedObserver)
    NamedTuple(prop => getproperty(O,prop) for O in Observer.Observers for prop in propertynames(O))
end