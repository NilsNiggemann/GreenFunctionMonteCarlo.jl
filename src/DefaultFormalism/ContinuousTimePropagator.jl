struct ContinuousTimePropagator{T<:AbstractFloat} <: AbstractPropagator 
    dτ::T
end
ContinuousTimePropagator(dτ::Real) = ContinuousTimePropagator(float(dτ))

function performMarkovStep!(x::AbstractConfig,moveWeights::AbstractVector,H::AbstractSignFreeOperator)
    moveidx = StatsBase.sample(StatsBase.Weights(moveWeights))
    move = get_move(H,moveidx)
    apply!(x,move)
    return moveidx
end
