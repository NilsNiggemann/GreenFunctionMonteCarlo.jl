"""
    AbstractGuidingFunction

An abstract type for all guiding function implementations in Green Function Monte Carlo. Subtypes of this correspond to implementations of concrete guiding functions, particularly variational wavefunctions. For optimal performance, it is recommended to implement the function log_psi_diff! for the specific guiding function type.
# Interface
- `logψ(x::AbstractArray)`: return the logarithm of the guiding function evaluated at the configuration `x`.
- `logψ(x::AbstractArray, H::AbstractHilbertSpace)`: return the logarithm of the guiding function evaluated at the configuration `x` in the specified `HilbertSpace`.
"""
abstract type AbstractGuidingFunction end
abstract type AbstractGuidingFunctionBuffer end

struct EmptyGWFBuffer <: AbstractGuidingFunctionBuffer end

@inline (logψ::AbstractGuidingFunction)(x::AbstractArray) = throw(MethodError(logψ, (x,)))
@inline (logψ::AbstractGuidingFunction)(x::AbstractArray,::AbstractHilbertSpace) = logψ(x)

function log_psi_diff(x::AbstractArray, dx::AbstractArray,logψ::AbstractGuidingFunction, Buffer::AbstractGuidingFunctionBuffer,Hilbert::AbstractHilbertSpace)
    logpsi = logψ(x)
    apply_move!(x,dx)
    if fulfills_constraint(x,Hilbert) 
        logpsi_new = logψ(x)
    else 
        logpsi_new = -Inf
    end
    apply_inverse!(x,dx)
    return logpsi_new - logpsi
end
@inline log_psi_diff(x::AbstractArray, dx::AbstractArray,logψ::AbstractGuidingFunction,Hilbert::AbstractHilbertSpace) = log_psi_diff(x,dx,logψ,EmptyGWFBuffer(),Hilbert)

function log_psi_diff!(logW::AbstractMatrix, X::AbstractWalkerEnsemble, dX::AbstractArray,logψ::AbstractGuidingFunction, Buffer::AbstractGuidingFunctionBuffer,Hilbert::AbstractHilbertSpace) 
    for j in axes(dX,2)
        dx = @view dX[:,j]
        for (α,x) in enumerate(eachConfig(X))
            logW[α,j] = log_psi_diff(x,dx,logψ,Buffer,Hilbert)
        end
    end
    return logW
end
@inline log_psi_diff!(logW::AbstractMatrix, X::AbstractWalkerEnsemble, dX::AbstractArray,logψ::AbstractGuidingFunction,Hilbert::AbstractHilbertSpace) = log_psi_diff!(logW,X,dX,logψ,EmptyGWFBuffer(),Hilbert)