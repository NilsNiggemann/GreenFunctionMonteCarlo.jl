struct KernelParallelization <: AbstractParallelizationScheme
    MemoryType::UnionAll
end
import AcceleratedKernels as AK


struct MatrixMoveOperator{MoveType,T,DiagType} <: AbstractSignFreeOperator
    moves::Matrix{MoveType}
    off_diag::Vector{T}
    diag::DiagType
end

@inline get_move(O::MatrixMoveOperator, idx::Integer) = @view O.moves[:,idx]
@inline get_diagonal(O::MatrixMoveOperator) = O.diag
@inline get_offdiagonal_elements(O::MatrixMoveOperator) = O.off_diag


struct WaveFunctionValues{T,Arr1<:AbstractArray{T,1},Arr2<:AbstractArray{T,2}}
    psix_alpha::Arr1
    psi_m_alpha::Arr2
end

struct ArrayWalkerEnsemble{T,ConfMat<:AbstractMatrix{T},ConnConfType<:AbstractArray{T,3},WeightMatType<:AbstractMatrix,BuffType<:WaveFunctionValues} <: AbstractWalkerEnsemble
    Configs::ConfMat
    ConnConfigs::ConnConfType
    WalkerWeights::Vector{Float64}
    MoveWeights::WeightMatType
    reconfigurationList::Vector{Int}
    local_energies::Vector{Float64}
    GWFBuffers::BuffType
end

getConfig(X::ArrayWalkerEnsemble, α) = @view X.Configs[:, α]
getConfigs(X::ArrayWalkerEnsemble) = X.Configs
getMoveWeights(X::ArrayWalkerEnsemble, α) = @view X.MoveWeights[:, α]
getMoveWeights(X::ArrayWalkerEnsemble) = X.MoveWeights
getWalkerWeights(X::ArrayWalkerEnsemble) = X.WalkerWeights
getBuffer(X::ArrayWalkerEnsemble) = X.GWFBuffers
getReconfigurationList(X::ArrayWalkerEnsemble) = X.reconfigurationList
getLocalEnergies(X::ArrayWalkerEnsemble) = X.local_energies

function allocate_walkerEnsemble(conf::AbstractVector{T}, logψ::ParametrizedFunction,NWalkers::Integer,numMoves::Integer,parallelization::KernelParallelization) where {T}
    MemoryType = parallelization.MemoryType
    # parentConf = parent(conf)

    Nsites = length(conf)
    configs = MemoryType{T}(undef, Nsites, NWalkers)
    AK.foraxes(configs,2) do α
        @. configs[:, α] = conf
    end

    conn_configs = MemoryType{T}(undef, Nsites, numMoves, NWalkers) .= 0

    weights = zeros(NWalkers)
    move_weights = MemoryType{Float64}(undef,numMoves, NWalkers) .= 0
    # reconfigurationList = MemoryType{Int}(NWalkers)
    reconfigurationList = zeros(Int, NWalkers)
    # local_energies = MemoryType{Float64}(undef,NWalkers) .= 0
    local_energies = zeros(NWalkers)

    psi_return = logψ(conf)

    psi_x_alpha = MemoryType{typeof(psi_return)}(undef, NWalkers) .= 0
    psi_m_alpha = MemoryType{typeof(psi_return)}(undef, numMoves, NWalkers) .= 0
    Buffers = WaveFunctionValues(psi_x_alpha, psi_m_alpha)

    return ArrayWalkerEnsemble(configs, conn_configs, weights, move_weights, reconfigurationList, local_energies, Buffers)
end

@views function getOffdiag!(ConnConfigs::AbstractArray{T,3}, Configs::AbstractMatrix{T}, moveMatrix::AbstractMatrix) where {T<:Integer}
    AK.foraxes(Configs,2) do α
        for m in axes(moveMatrix,2)
            @. ConnConfigs[:, m, α] = Configs[:, α] + moveMatrix[:, m]
        end
    end
    return ConnConfigs
end
function getOffdiag!(ConnConfigs::AbstractArray{Bool,3}, Configs::AbstractMatrix{Bool}, moveMatrix::AbstractMatrix{Bool})
    AK.foraxes(Configs,2) do α
        for m in axes(moveMatrix,2)
            for i in axes(Configs,1)
                ConnConfigs[i, m, α] = xor(Configs[i, α], moveMatrix[i, m])
            end
        end
    end
    return ConnConfigs
end

"""
    map_function!(out_i, xij, func, params)

Applies the `func` function to each element in `xij` with the given `params`, storing the results in `out_i`.

# Arguments
- `out_i`: Output array to store the results (modified in-place).
- `xij`: Input data, typically an array or matrix of coordinates or values.
- `func`: Function to be mapped over `xij`, usually representing a wavefunction or similar operation.
- `params`: Additional parameters required by `func`.

# Notes
This function is intended for use in GPU-based computations, and may leverage GPU acceleration for performance.
"""
function map_function!(out_i, xij, func, params)
    AK.foreachindex(out_i) do m_alpha
        out_i[m_alpha] = func(view(xij, :, m_alpha), params)
    end
    return out_i
end
map_function!(out_i, xij, func::ParametrizedFunction) = map_function!(out_i, xij, func.logpsi, get_params(func))
map_function!(out_i, xij, func) = error("map_function! not implemented for function type $(typeof(func)). Use ParametrizedFunction or map_function!(out_i, xij, func, params).")

function get_psi_x!(WE::ArrayWalkerEnsemble, logψ::ParametrizedFunction)
    Configs = getConfigs(WE)
    Buffer = getBuffer(WE)
    psi_x = map_function!(Buffer.psix_alpha, Configs, logψ)
    return psi_x
end

function get_psi_m_alpha!(WE::ArrayWalkerEnsemble, logψ::ParametrizedFunction)
    ConnConfigs = WE.ConnConfigs
    Buffer = getBuffer(WE)
    Nsites, numMoves, Nwalkers = size(ConnConfigs)
    psi_m_alpha = Buffer.psi_m_alpha
    psi_m_alpha_reshape = reshape(psi_m_alpha, numMoves*Nwalkers)
    ConnConfigs_reshape = reshape(ConnConfigs, Nsites, numMoves*Nwalkers)

    map_function!(psi_m_alpha_reshape, ConnConfigs_reshape, logψ)
    return psi_m_alpha
end
getOffdiag!(WE::ArrayWalkerEnsemble, moveMatrix::AbstractMatrix) = getOffdiag!(WE.ConnConfigs, WE.Configs, moveMatrix)


function _get_markov_weights!(weights::AbstractMatrix,psi_x,psi_m_alpha,H::MatrixMoveOperator,Hilbert::AbstractHilbertSpace)
    Hxy = get_offdiagonal_elements(H)
    AK.foraxes(weights,2) do alpha
        for m in axes(weights,1)
            weights[m, alpha] = -Hxy[m]*exp(psi_m_alpha[m, alpha] - psi_x[alpha])
        end
    end
    return weights
end
function get_markov_weights!(WE::ArrayWalkerEnsemble,H::MatrixMoveOperator,logψ::ParametrizedFunction, Hilbert)
    getOffdiag!(WE, H.moves)
    psi_x = get_psi_x!(WE, logψ)
    psi_m_alpha = get_psi_m_alpha!(WE, logψ)
    return _get_markov_weights!(WE.MoveWeights, psi_x, psi_m_alpha, H, Hilbert)
end

function continuos_time_propagation!(WE::ArrayWalkerEnsemble, H::AbstractSignFreeOperator, logψ, Hilbert::AbstractHilbertSpace, dτ::Real, w_avg_estimate::Real, parallelization::KernelParallelization, RNG::Random.AbstractRNG = Random.default_rng())
    Configs = getConfigs(WE)
    ConnConfigs = getOffdiag!(WE, H.moves)
    
    psi_x = map_function!(similar(Configs), Configs, logψ, nothing)

    moveWeights = getMoveWeights(WE)
    Nwalkers = size(Configs, 2)
    # getMoveWeights!(moveWeights, WE, H, params)
    log_w = zeros(Nwalkers)
    H_xx = zeros(Nwalkers)
    
    Hxx = get_diagonal(H)
    map_function!(H_xx, Configs, Hxx)

    get_markov_weights!(WE, H, logψ, Hilbert)

    for α in eachindex(WE)
        log_w = 0.0
        get_markov_weights!(moveWeights, Config, H, logψ, Hilbert, GWFBuffer)
        Hxx = get_diagonal(H)
        H_xx = Hxx(Config)
        el_x = H_xx + getLocalEnergy(moveWeights)
        τleft = dτ
        while τleft > 0
            ξ = rand(RNG)
            τ_step = min(τleft, log(1 - ξ) / (el_x - H_xx))
            τleft -= τ_step
            if isinf(τleft)
                @info "" τ_step el_x H_xx τleft maximum(moveWeights)
                error("Infinite propagation time encountered. Check for too large values in guiding wavefunction or its Buffers!")
            end
            log_w += -τ_step * el_x
            if τleft > 0 
                last_move = performMarkovStep!(Config, moveWeights, H, RNG)
                post_move_affect!(GWFBuffer, Config, last_move, logψ)
                get_markov_weights!(moveWeights, Config, H, logψ, Hilbert, GWFBuffer)

                H_xx = Hxx(Config)
                el_x = H_xx + getLocalEnergy(moveWeights)
            end
        end
        w = exp(log_w - dτ * w_avg_estimate)
        WalkerWeights = getWalkerWeights(WE)
        WalkerWeights[α] = w
        localEnergies = getLocalEnergies(WE)
        localEnergies[α] = el_x
    end
    return nothing
end

