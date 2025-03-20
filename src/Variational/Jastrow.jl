"""
    struct Jastrow{T<:Real} <: AbstractGuidingFunction

A structure representing the Jastrow factor in a variational Monte Carlo simulation. This function strikes a good balance between accuracy and computational efficiency in the context of Green Function Monte Carlo.
``
log(ψ(x)) = \\sum_i m_i x_i + \\frac{1}{2}\\sum_{i,j} v_{ij} x_i x_j
``

# Type Parameters
- `T<:Real`: The numeric type used for the parameters of the Jastrow factor 
  (e.g., `Float64`, `Float32`).

# Usage
logψ = Jastrow(N,Float64) # Create a Jastrow factor for `N` sites with parameters of type `Float64`.
logψ = Jastrow(conf,Float64) # Create a Jastrow factor for configurations similar to `conf`.
"""
struct Jastrow{T<:Real} <: AbstractGuidingFunction
    m_i::Vector{T}
    v_ij::Matrix{T}
    buffer_reset_max::Int
end

function Jastrow(N::Int,Type = Float32; buffer_reset_max=1000)
    m_i = zeros(Type,N)
    v_ij = zeros(Type,N,N)
    return Jastrow(m_i,v_ij,buffer_reset_max)
end

Jastrow(conf::AbstractArray,args...;kwargs...) = Jastrow(length(conf),args...;kwargs...)

get_m_i(logψ::Jastrow) = logψ.m_i
get_v_ij(logψ::Jastrow) = logψ.v_ij

function get_params(logψ::Jastrow)
    return RecursiveArrayTools.ArrayPartition(logψ.m_i,get_v_ij(logψ))
end

function (logψ::Jastrow)(x::AbstractConfig)
    m = get_m_i(logψ)
    v = get_v_ij(logψ)
    evaluate_jastrow(x,m,v)
end

guidingfunc_name(F::Jastrow) = "Jastrow"

function evaluate_jastrow(x,m::AbstractVector,v::AbstractMatrix)
    log_m = zero(eltype(m))
    for i in eachindex(IndexLinear(),x)
        log_m += m[i] * x[i]
    end
    x_lin = reshape(x, length(x))
    log_v = 0.5*LinearAlgebra.dot(x_lin,v,x_lin)

    return log_m + log_v
end

struct SimpleJastrow_GWF_Buffer{T<:Number}
    x_i::Vector{T}
    h_i::Vector{T}
    buffer_reset_counter::Ref{Int}
    # prefac_moves::Matrix{T}
end

function setBuffer!(A::SimpleJastrow_GWF_Buffer,B::SimpleJastrow_GWF_Buffer)
    A.x_i .= B.x_i
    A.h_i .= B.h_i
    A.buffer_reset_counter[] = B.buffer_reset_counter[]
    # A.prefac_moves .= B.prefac_moves
    return A
end


function allocate_GWF_buffer(logψ::Jastrow{T},conf) where T
    h_i = zeros(T,length(conf))
    buffer_reset_counter = Ref(0)
    # prefac_moves = _precompute_prefac_moves(logψ,S)
    x_i = zeros(T,length(conf))
    Buff =  SimpleJastrow_GWF_Buffer(x_i,h_i,buffer_reset_counter)
    compute_GWF_buffer!(Buff,logψ,conf)
    return Buff
end
function pre_move_affect!(Buffer::SimpleJastrow_GWF_Buffer,conf::AbstractConfig,logψ::Jastrow) 
    return Buffer
end

Base.@propagate_inbounds function compute_GWF_buffer!(Buffer::SimpleJastrow_GWF_Buffer,logψ::Jastrow,x::AbstractConfig)

    v = get_v_ij(logψ)
    h = Buffer.h_i
    m = get_m_i(logψ)
    @boundscheck checkbounds(x,axes(v,1))
    @boundscheck checkbounds(h,axes(v,1))

    # LoopVectorization.@turbo for i in axes(v,1)
    for i in axes(v,1)
        hi = m[i]
        for j in axes(v,2)
            hi += x[j]*v[j,i]
        end
        h[i] = hi
    end
    return Buffer
end


function post_move_affect!(Buffer::SimpleJastrow_GWF_Buffer,logψ::Jastrow,Config::AbstractConfig,move::AbstractMove)
    sites = affected_sites(move)
    dx = move_dx(move,Config)

    v = get_v_ij(logψ)
    h = Buffer.h_i
    
    # LoopVectorization.@turbo for (index,i) in enumerate(inds)
    buffer_reset_counter = Buffer.buffer_reset_counter
    buffer_reset_max = logψ.buffer_reset_max
    if buffer_reset_counter[] >= buffer_reset_max
        buffer_reset_counter[] = 0
        compute_GWF_buffer!(Buffer,logψ,Config)
        return Buffer
    end

    for (idx,i) in enumerate(sites)
        s = dx[idx]
        
        # LoopVectorization.@turbo for j in eachindex(h)
        for j in eachindex(h)
            h[j] += v[j,i]*s
        end
    end

    return Buffer
end

@inline function log_psi_diff(Config::AbstractConfig, move::AbstractMove, logψ::Jastrow{T}, Buffer::SimpleJastrow_GWF_Buffer, Hilbert::AbstractHilbertSpace) where T
    isapplicable(Config,move, Hilbert) || return -Inf

    (;h_i) = Buffer
    
    sites = affected_sites(move)
    dx = move_dx(move,Config)

    log_h = zero(T)
    # m = get_m_i(logψ)
    # @inbounds @simd for idx in affected_sites
    for (idx,i) in enumerate(sites)
        s = dx[idx]
        log_h += h_i[i]*s #note: h_i contains m_i
    end

    vij = get_v_ij(logψ)
    for (i_k,k) in enumerate(sites)
        for (i_j,j) in enumerate(sites)
            log_h += 0.5*dx[i_k]*dx[i_j]*vij[j,k]
        end
    end
    return log_h
end