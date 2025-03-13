struct ContinuousTimePropagator{T<:AbstractFloat} <: AbstractPropagator 
    dτ::T
end
ContinuousTimePropagator(dτ::Real) = ContinuousTimePropagator(float(dτ))
