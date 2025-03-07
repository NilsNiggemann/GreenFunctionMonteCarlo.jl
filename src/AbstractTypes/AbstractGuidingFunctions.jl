abstract type AbstractGuidingFunction end
abstract type AbstractGuidingFunctionBuffer end

(logψ::AbstractGuidingFunction)(x::AbstractConfig) = throw(MethodError(logψ, (x,)))
(logψ::AbstractGuidingFunction)(x::AbstractConfig,::AbstractHilbertSpace) = logψ(x)

log_psi_diff(logψ::AbstractGuidingFunction, x::AbstractConfig, y::AbstractConfig, Buffer::AbstractGuidingFunctionBuffer) = logψ(x) - logψ(y)