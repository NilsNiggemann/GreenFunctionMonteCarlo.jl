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
Jastrow(m_i::AbstractVector,v_ij::AbstractMatrix; buffer_reset_max=3_000) = Jastrow(m_i,_make_symmetric(v_ij),buffer_reset_max)
Jastrow(conf::AbstractArray,args...;kwargs...) = Jastrow(length(conf),args...;kwargs...)

_make_symmetric(v::AbstractMatrix) = copy(v) .= LinearAlgebra.Symmetric(v)

function Jastrow(N::Int,Type = Float32; buffer_reset_max=3_000)
    m_i = zeros(Type,N)
    v_ij = zeros(Type,N,N)
    return Jastrow(m_i,v_ij,buffer_reset_max)
end


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

struct SimpleJastrow_GWF_Buffer{T<:Number} <: AbstractGuidingFunctionBuffer
    x_i::Vector{T}
    h_i::Vector{T}
    buffer_reset_counter::Array{Int,0}
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
    buffer_reset_counter = zeros(Int)
    # prefac_moves = _precompute_prefac_moves(logψ,S)
    x_i = zeros(T,length(conf))
    Buff =  SimpleJastrow_GWF_Buffer(x_i,h_i,buffer_reset_counter)
    compute_GWF_buffer!(Buff,logψ,conf)
    return Buff
end
function pre_move_affect!(Buffer::SimpleJastrow_GWF_Buffer,conf::AbstractConfig,logψ::Jastrow) 
    return Buffer
end

Base.@propagate_inbounds function compute_GWF_buffer!(Buffer::SimpleJastrow_GWF_Buffer,logψ::Jastrow,config::AbstractConfig)

    v = get_v_ij(logψ)
    h = Buffer.h_i
    m = get_m_i(logψ)
    x = parent(config)
    @boundscheck checkbounds(x,axes(v,1))
    @boundscheck checkbounds(h,axes(v,1))
    @boundscheck checkbounds(h,axes(v,2))

    LoopVectorization.@turbo for i in axes(v,1)
        hi = m[i]
        for j in axes(v,2)
            hi += x[j]*v[j,i]
        end
        h[i] = hi
    end
    return Buffer
end


function post_move_affect!(Buffer::SimpleJastrow_GWF_Buffer,Config::AbstractConfig,move::AbstractMove,logψ::Jastrow)
    sites = affected_sites(move)
    dx = move_dx_after(move,Config)

    v = get_v_ij(logψ)
    h = Buffer.h_i
    buffer_reset_counter = Buffer.buffer_reset_counter
    buffer_reset_max = logψ.buffer_reset_max
    if buffer_reset_counter[] >= buffer_reset_max # recompute buffers only occasionally to avoid accumulation of floating point errors 
        buffer_reset_counter[] = 0
        compute_GWF_buffer!(Buffer,logψ,Config)
        return Buffer
    end
    
    for (idx,i) in enumerate(sites)
        s = dx[idx]
        
        LoopVectorization.@turbo for j in eachindex(h)
        # for j in eachindex(h)
            h[j] += v[j,i]*s
        end
    end
    buffer_reset_counter[] += 1
    return Buffer
end

@inline function log_psi_diff(Config::AbstractConfig, move::AbstractMove, logψ::Jastrow{T}, Buffer::SimpleJastrow_GWF_Buffer, Hilbert::AbstractHilbertSpace) where T
    isapplicable(Config,move, Hilbert) || return -Inf
    (;h_i) = Buffer
    
    sites = affected_sites(move)
    dx = move_dx_jastrow(move,Config)

    vij = get_v_ij(logψ)

    return _jastrow_diff_kernel(sites,dx,h_i,vij)
end
# hack into move_dx to allow even faster path which utilizes booleans
@inline function move_dx_jastrow(move::FlipMove{<:SA.SVector},x::AbstractConfig{Bool})
    map(move.inds) do i
        x[i]
    end
end
# default fallback, use move_dx which returns a Number type
@inline function move_dx_jastrow(move::AbstractMove,x::AbstractConfig)
    move_dx_before(move,x)
end

const LV_COMPATIBLE_VECTOR = Union{Vector,SA.SVector}
const LV_COMPATIBLE_MATRIX = Union{Matrix,SA.SMatrix}

# fast compile-time branch for types allowed by LoopVectorization
@inline function _jastrow_diff_kernel(sites::LV_COMPATIBLE_VECTOR,dx::SA.SVector{N,Bool},h_i::LV_COMPATIBLE_VECTOR,vij::LV_COMPATIBLE_MATRIX) where N
    log_h = zero(eltype(h_i))
    # note: dx[idx] is a boolean. false corresponds to x=0 => dx = 1, true corresponds to x=1 => dx = -1
    LoopVectorization.@turbo for idx in eachindex(sites)
        i = sites[idx]
        s = dx[idx]
        his = ifelse(s,-h_i[i],h_i[i])
        log_h += his
        # log_h += h_i[i]*s #note: h_i contains m_i
    end

    log_h2 = zero(log_h)
    LoopVectorization.@turbo for (i_k,k) in enumerate(sites) #todo: precompute this for all moves in the GWF Buffer
        s1 = dx[i_k]
        for (i_j,j) in enumerate(sites)
            # log_h += 0.5*dx[i_k]*dx[i_j]*vij[j,k]
            vv = vij[j,k]
            s2 = dx[i_j]
            ssv = ifelse(s1 ⊻ s2,-vv,vv)
            log_h2 += ssv
        end
    end
    return log_h + 0.5log_h2
end

@inline function _jastrow_diff_kernel(sites::LV_COMPATIBLE_VECTOR,dx::LV_COMPATIBLE_VECTOR,h_i::LV_COMPATIBLE_VECTOR,vij::LV_COMPATIBLE_MATRIX)
    log_h = zero(eltype(h_i))
    LoopVectorization.@turbo for idx in eachindex(sites)
        i = sites[idx]
        s = dx[idx]
        log_h += h_i[i]*s #note: h_i contains m_i
    end

    log_h2 = zero(eltype(h_i))
    LoopVectorization.@turbo for (i_k,k) in enumerate(sites) #todo: precompute this for all moves in the GWF Buffer
        s1 = dx[i_k]
        sumterm = zero(log_h2)
        for (i_j,j) in enumerate(sites)
            # log_h2 += dx[i_k]*dx[i_j]*vij[j,k]
            sumterm += dx[i_j]*vij[j,k]
        end
        log_h2 += s1*sumterm
    end
    return log_h + 0.5log_h2
end
# slower, generic fallback
@inline function _jastrow_diff_kernel(sites,dx::AbstractVector{<:Number},h_i,vij)
    log_h = zero(eltype(h_i))

    @inbounds @simd for idx in eachindex(sites)
        i = sites[idx]
        s = dx[idx]
        log_h += h_i[i]*s #note: h_i contains m_i
    end

    log_h2 = zero(log_h)
    @inbounds @simd for i_k in eachindex(sites)
        k = sites[i_k]
        s1 = dx[i_k]
        sumterm = zero(log_h)
        for i_j in eachindex(sites)
            j = sites[i_j]
            # log_h += 0.5*dx[i_k]*dx[i_j]*vij[j,k]
            sumterm += dx[i_j]*vij[j,k]
        end
        log_h2 += s1*sumterm
    end
    return log_h + 0.5log_h2
end