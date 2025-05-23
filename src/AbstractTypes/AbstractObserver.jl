"""
Abstract supertype for observables in GFMC which are recorded during the run. Observables should always contain the table recording reconfiguration processes, the energies, and the total weights of each Markov step.
# Interface: 
- `saveObservables_before!(Observables,i,Walkers,propagator)`: Saves the observables for the given iteration `i` and walker ensemble `Walkers`.
- `saveObservables_after!(Observables,i,Walkers,propagator)`: Saves the observables for the given iteration `i` and walker ensemble `Walkers` after reconfiguration.
"""
abstract type AbstractObserver end

saveObservables_before!(Observables::AbstractObserver,i,Walkers,H,reconfiguration) = nothing

saveObservables_after!(Observables::AbstractObserver,i,Walkers,H,reconfiguration) = nothing

"""
Abstract supertype for observables which are diagonal in the computational basis and may be measured for free in GFMC. An observable must be a function that takes a configuration and returns an array.
# Interface:
- `obs(::O)` returns the buffer array for the output of the observable.
- `O(out,Conf)`: Writes the observable for the given configuration to the preallocated array `out`. Returns out.
- `Base.copy(O::O)`: Returns a copy of the observable.

If no preallocated array is given, the observable defaults to using the buffer array `obs(O)`.
"""
abstract type AbstractDiagObservable end