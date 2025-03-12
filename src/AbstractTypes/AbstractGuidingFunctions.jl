"""
    AbstractGuidingFunction

An abstract type for all guiding function implementations in Green Function Monte Carlo. Subtypes of this correspond to implementations of concrete guiding functions, particularly variational wavefunctions. For optimal performance, it is recommended to implement the function log_psi_diff! for the specific guiding function type.
# Interface
- `logψ(x::AbstractArray)`: return the logarithm of the guiding function evaluated at the configuration `x`.
- `logψ(x::AbstractArray, H::AbstractHilbertSpace)`: return the logarithm of the guiding function evaluated at the configuration `x` in the specified `HilbertSpace`.
- `log_psi_diff(x::AbstractArray, dx::AbstractArray, logψ::AbstractGuidingFunction, Buffer::AbstractGuidingFunctionBuffer, Hilbert::AbstractHilbertSpace)`: return the logarithm of the ratio of the guiding function evaluated at the configuration `x` and `x+dx` in the specified `HilbertSpace`. Returns `-Inf` if the move is not applicable.
- get_params(logψ::AbstractGuidingFunction): return the parameters of the guiding function as a linear Array. It is recommended to use RecursiveArrayTools.jl for this purpose.
- HDF5.h5write(file::AbstractString, name::AbstractString, logψ::AbstractGuidingFunction): write the guiding function to an HDF5 file.
- allocate_GWF_buffers(logψ::AbstractGuidingFunction, NBuffers::Integer): allocate NBuffers instances of a buffer for the guiding function. Defaults to an array of EmptyGWFBuffer instances.
- compute_GWF_buffer!(Buffer::AbstractGuidingFunctionBuffer, logψ::AbstractGuidingFunction, x): compute the full buffer for the guiding function at the configuration `x`. 
- pre_move_affect!(Buffer::AbstractGuidingFunctionBuffer, x, dx, logψ::AbstractGuidingFunction): perform any necessary operations before the move is applied.
- post_move_affect!(Buffer::AbstractGuidingFunctionBuffer, x, dx, logψ::AbstractGuidingFunction): perform any necessary operations after the move is applied.
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

get_params(logψ::AbstractGuidingFunction) = throw(MethodError(get_params, (logψ,)))
HDF5.h5write(file::AbstractString,name::AbstractString,logψ::AbstractGuidingFunction) = throw(MethodError(h5write, (file,name,logψ)))

allocate_GWF_buffers(logψ::AbstractGuidingFunction, NBuffers::Integer) = fill(EmptyGWFBuffer(),NBuffers)

compute_GWF_buffer!(Buffer::AbstractGuidingFunctionBuffer,logψ::AbstractGuidingFunction,x) = nothing

function compute_GWF_buffers!(Walkers::AbstractWalkerEnsemble,logψ::AbstractGuidingFunction)
    Buffers = getBuffers(Walkers)
    X = getConfigs(Walkers)
    Threads.@threads for i in eachindex(Buffers,X)
        compute_GWF_buffer!(Buffers[i],X[i],logψ)
    end
    return nothing
end

function pre_move_affect!(Buffer::AbstractGuidingFunctionBuffer,x,dx,logψ::AbstractGuidingFunction)
    return nothing
end

function post_move_affect!(Buffer::AbstractGuidingFunctionBuffer,x,dx,logψ::AbstractGuidingFunction)
    return nothing
end