"""
Abstract supertype for observables. An observable must be a function that takes a spin configuration and returns an array.
To define a subtype of O of `AbstractObservable`, one must define the following functions:

- obs(::O) returns the buffer array for the output of the observable.
- `O(out,Conf)`: Writes the observable for the given configuration to the preallocated array `out`. Returns out.
- Base.copy(O::O): Returns a copy of the observable.

If no preallocated array is given, the observable defaults to using the buffer array `obs(O)`.
"""
abstract type AbstractObservable end
obs_size(O::AbstractObservable) = size(obs(O))
(O::AbstractObservable)(Conf) = O(obs(O),Conf)
Base.show(io::IO,::MIME"text/plain",O::Obs) where {Obs <: AbstractObservable} = print(io, "$Obs,", " ∈ ", obs_size(O))
Base.display(io::IO,::MIME"text/plain",O::Obs) where {Obs <: AbstractObservable} = print(io, "$Obs,", " ∈ ", obs_size(O))