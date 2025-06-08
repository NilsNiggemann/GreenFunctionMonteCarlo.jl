"""
    AbstractGuidingFunction

An abstract type for all guiding function implementations in Green Function Monte Carlo. Subtypes of this correspond to implementations of concrete guiding functions, particularly variational wavefunctions. For optimal performance, it is recommended to implement the function log_psi_diff! for the specific guiding function type.
# Interface
- `logψ(x::AbstractArray)`: return the logarithm of the guiding function evaluated at the configuration `x`.
- `logψ(x::AbstractArray, H::AbstractHilbertSpace)`: return the logarithm of the guiding function evaluated at the configuration `x` in the specified `HilbertSpace`.
- `log_psi_diff(x::AbstractArray, dx::AbstractArray, logψ::AbstractGuidingFunction, Buffer::AbstractGuidingFunctionBuffer, Hilbert::AbstractHilbertSpace)`: return the logarithm of the ratio of the guiding function evaluated at the configuration `x` and `x+dx` in the specified `HilbertSpace`. Returns `-Inf` if the move is not applicable.
- `get_params(logψ::AbstractGuidingFunction)`: return the parameters of the guiding function as a linear Array. It is recommended to use RecursiveArrayTools.jl for this purpose.
# Interface (optional)
- `allocate_GWF_buffers(logψ::AbstractGuidingFunction, NBuffers::Integer)`: allocate NBuffers instances of a buffer for the guiding function. Defaults to an array of EmptyGWFBuffer instances.
- `compute_GWF_buffer!(Buffer::AbstractGuidingFunctionBuffer, logψ::AbstractGuidingFunction, x)`: compute the full buffer for the guiding function at the configuration `x`. 
- `pre_move_affect!(Buffer::AbstractGuidingFunctionBuffer, x, moves, logψ::AbstractGuidingFunction)`: perform any necessary operations before computing ratios ψx´_ψx.
- `post_move_affect!(Buffer::AbstractGuidingFunctionBuffer, x, dx, logψ::AbstractGuidingFunction)`: perform any necessary operations after the move is applied.
- `HDF5.h5write(file::AbstractString, name::AbstractString, logψ::AbstractGuidingFunction)`: write the guiding function to an HDF5 file.
"""
abstract type AbstractGuidingFunction end

"""
    AbstractGuidingFunctionBuffer

An abstract type that serves as a base for defining guiding function buffer structures. 

# Interface
- `setBuffer!(BA::AbstractGuidingFunctionBuffer, BB::AbstractGuidingFunctionBuffer)`: set the buffer `BA` to the values of the buffer `BB`.
"""
abstract type AbstractGuidingFunctionBuffer end

struct EmptyGWFBuffer <: AbstractGuidingFunctionBuffer end
struct NotImplementedBuffer <: AbstractGuidingFunctionBuffer end

function (logψ::AbstractGuidingFunction) end
@inline (logψ::AbstractGuidingFunction)(x::AbstractArray,::AbstractHilbertSpace) = logψ(x)

function log_psi_diff(x, dx::AbstractMove, logψ::AbstractGuidingFunction,logψx::AbstractFloat,Hilbert::AbstractHilbertSpace)
    apply!(x,dx)
    if fulfills_constraint(x,Hilbert) 
        logpsi_new = logψ(x)
    else 
        logpsi_new = -Inf
    end
    apply_inverse!(x,dx)
    return logpsi_new - logψx
end

setBuffer!(BA::EmptyGWFBuffer,BB::EmptyGWFBuffer) = BB
setBuffer!(BA::NotImplementedBuffer,BB::NotImplementedBuffer) = BB

function get_params end
# HDF5.h5write(file::AbstractString,name::AbstractString,logψ::AbstractGuidingFunction) = throw(MethodError(h5write, (file,name,logψ)))

"""
    allocate_GWF_buffers(logψ::AbstractGuidingFunction, NBuffers::Integer, x)

Allocates a specified number of guiding wave function (GWF) buffers.

# Arguments
- `logψ::AbstractGuidingFunction`: The guiding function for which the buffers are being allocated.
- `NBuffers::Integer`: The number of buffers to allocate.
- `x`: An exemplary configuration

# Returns
- An array filled with `NBuffers` instances of `AbstractGuidingFunctionBuffer`. Defaults to an array of `NotImplementedBuffer` instances.
"""
function allocate_GWF_buffers(logψ::AbstractGuidingFunction, x, NBuffers)
    Buff1 = allocate_GWF_buffer(logψ,x)
    allBuffs = Vector{typeof(Buff1)}(undef,NBuffers)
    Threads.@threads for i in 1:NBuffers
        allBuffs[i] = allocate_GWF_buffer(logψ,x)
    end
    return allBuffs
end


compute_GWF_buffer!(Buffer::AbstractGuidingFunctionBuffer,logψ::AbstractGuidingFunction,x) = Buffer

function compute_GWF_buffers!(Walkers::AbstractWalkerEnsemble,logψ::AbstractGuidingFunction)
    Threads.@threads for α in eachindex(Walkers)
        Buffer = getBuffer(Walkers,α)
        x = getConfig(Walkers,α)
        compute_GWF_buffer!(Buffer,logψ,x)
    end
end

"""
    Updates the guiding function buffer before any move is applied.
"""
function pre_move_affect!(Buffer::AbstractGuidingFunctionBuffer,x,logψ::AbstractGuidingFunction)
    return Buffer
end

function post_move_affect!(Buffer::AbstractGuidingFunctionBuffer,x,dx,logψ::AbstractGuidingFunction)
    return Buffer
end