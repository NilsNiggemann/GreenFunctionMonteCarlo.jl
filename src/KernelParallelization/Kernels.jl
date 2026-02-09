struct KernelParallelization{MemoryType} <: AbstractParallelizationScheme
end
KernelParallelization(::Type{T}) where {T<:AbstractArray} = KernelParallelization{T}()
getMemoryType(par::KernelParallelization{T}) where T = T

struct MatrixMoveOperator{MatType,VecType,DiagType} <: AbstractSignFreeOperator
    moves::MatType
    off_diag::VecType
    diag::DiagType
end
function Adapt.adapt_structure(to, x::MatrixMoveOperator)
    diag = Adapt.adapt(to, x.diag)
    off_diag = Adapt.adapt(to, x.off_diag)
    moves = Adapt.adapt(to, x.moves)
    return MatrixMoveOperator(moves, off_diag, diag)
end

struct ArrayMove{T,Vectype<:AbstractVector{T}} <: AbstractMove
    move::Vectype
end
affected_sites(M::ArrayMove) = findall(!iszero, M.move)
apply!(x, M::ArrayMove) = (x .+= M.move)
apply!(x, M::ArrayMove{Bool}) = map!(xor, x, M.move)
move_dx(move::ArrayMove) = move.move

@inline get_move(O::MatrixMoveOperator, idx::Integer) = ArrayMove(@view O.moves[:,idx])
@inline get_diagonal(O::MatrixMoveOperator) = O.diag
@inline get_offdiagonal_elements(O::MatrixMoveOperator) = O.off_diag


struct WaveFunctionValues{T,Arr1<:AbstractArray{T,1},Arr2<:AbstractArray{T,2}}
    psix_alpha::Arr1
    psi_m_alpha::Arr2
end

struct ArrayWalkerEnsemble{T,ConfMat<:AbstractMatrix{T},ConnConfType<:AbstractArray{T,3},VecType<:AbstractVector,WeightMatType<:AbstractMatrix,BuffType<:WaveFunctionValues} <: AbstractWalkerEnsemble
    Configs::ConfMat
    ConnConfigs::ConnConfType
    WalkerWeights::VecType
    MoveWeights::WeightMatType
    reconfigurationList::Vector{Int}
    local_energies::VecType
    GWFBuffers::BuffType
end
Base.eachindex(X::ArrayWalkerEnsemble) = eachindex(X.WalkerWeights)
getConfig(X::ArrayWalkerEnsemble, α) = @view X.Configs[:, α]
getConfigs(X::ArrayWalkerEnsemble) = X.Configs
getMoveWeights(X::ArrayWalkerEnsemble, α) = @view X.MoveWeights[:, α]
getMoveWeights(X::ArrayWalkerEnsemble) = X.MoveWeights
getWalkerWeights(X::ArrayWalkerEnsemble) = X.WalkerWeights
getBuffer(X::ArrayWalkerEnsemble) = X.GWFBuffers
getBuffer(X::ArrayWalkerEnsemble,α) = NotImplementedBuffer()
getReconfigurationList(X::ArrayWalkerEnsemble) = X.reconfigurationList
getLocalEnergies(X::ArrayWalkerEnsemble) = X.local_energies

function set_to_config!(Configs, conf::AbstractVector{T}) where {T}
    AK.foraxes(Configs,2) do α
        Configs[:, α] .= conf
    end
    return Configs
end

function get_return_example(logψ::ParametrizedFunction, conf::AbstractVector)
    return AK.@allowscalar logψ(conf)
end

function allocate_walkerEnsemble(conf::AbstractVector{T}, logψ::ParametrizedFunction,NWalkers::Integer,numMoves::Integer,parallelization::KernelParallelization) where {T}
    MemoryType = getMemoryType(parallelization)
    # parentConf = parent(conf)
    conf_converted = MemoryType{T}(parent(conf))
    Nsites = length(conf)
    configs = MemoryType{T}(undef, Nsites, NWalkers)
    set_to_config!(configs, conf_converted)
    conn_configs = MemoryType{T}(undef, Nsites, numMoves, NWalkers) .= 0

    # weights = MemoryType{Float64}(undef,NWalkers) .= 0
    weights = zeros(NWalkers)
    move_weights = MemoryType{Float64}(undef,numMoves, NWalkers) .= 0
    # reconfigurationList = MemoryType{Int}(NWalkers)
    reconfigurationList = zeros(Int, NWalkers)
    # local_energies = MemoryType{Float64}(undef,NWalkers) .= 0
    # local_energies =  MemoryType{Float64}(undef,NWalkers) .= 0
    local_energies =  zeros(NWalkers)

    psi_return_type = typeof(get_return_example(logψ, conf_converted))

    isconcretetype(psi_return_type) || @warn ("Wavefunction return type is not a concrete type: $psi_return_type")
    psi_x_alpha = MemoryType{psi_return_type}(undef, NWalkers) .= 0
    psi_m_alpha = MemoryType{psi_return_type}(undef, numMoves, NWalkers) .= 0
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

function e_local_alpha!(e_local_alpha,Hxx_alpha,moveWeights)
    # Nwalkers = size(moveWeights, 2)
    # e_local_alpha .= Hxx_alpha
    # for α in eachindex(e_local_alpha)
    #     e_local_alpha[α] += getLocalEnergy(view(moveWeights, :, α))
    # end
    AK.foraxes(e_local_alpha) do alpha
        e_local_alpha[alpha] = Hxx_alpha[alpha] + getLocalEnergy(view(moveWeights, :, alpha))
    end
    return e_local_alpha
end

function get_move_idxs(moveWeights::Matrix,rng::Random.AbstractRNG)
    moveidxs = StatsBase.sample.(Ref(rng),StatsBase.Weights.(eachcol(moveWeights)))
    return moveidxs
end

function get_move_idxs(moveWeights::AbstractMatrix,rng::Random.AbstractRNG)
    # needed for now as StatsBase.sample does not work with GPU arrays
    Mat_reshape = reshape(moveWeights,:)
    Mat_matrix = Matrix(moveWeights)
    
    move_idxs = get_move_idxs(Mat_matrix,rng)
    ref_array = similar(moveWeights,Int, 0)

    Adapt.adapt(AK.get_backend(moveWeights), move_idxs)
end

function TEST_get_move_idxs(moveIdxs,moveWeights::AbstractMatrix,rng::Random.AbstractRNG)
    # needed for now as StatsBase.sample does not work with GPU arrays
    # Random.rand!(1:100,moveIdxs)
    moveIdxs .= 8
    return moveIdxs
end

# function performMarkovSteps!(Configs::AbstractMatrix,ConnectedConfs::AbstractArray{T,3},moveWeights::AbstractMatrix,walker_still_running,rng::Random.AbstractRNG) where {T}
#     # AK.@allowscalar moveidxs = StatsBase.sample.(Ref(rng),StatsBase.Weights.(eachcol(moveWeights)))
    
#     moveidxs = get_move_idxs(moveWeights,rng)
#     for α in eachindex(moveidxs)
#         if walker_still_running[α]
#             move = moveidxs[α]
#             Configs[:, α] .= view(ConnectedConfs, :, move, α)
#         end
#     end
#     return nothing
# end
function performMarkovSteps!(Configs::AbstractMatrix,ConnectedConfs::AbstractArray{T,3},moveWeights::AbstractMatrix,walker_still_running,rng::Random.AbstractRNG) where {T}
    moveidxs = get_move_idxs(moveWeights, rng)

    AK.foraxes(Configs,2) do α
        if walker_still_running[α]
            move = moveidxs[α]
            Configs[:, α] .= view(ConnectedConfs,:, move, α)
        end
    end
    return nothing
end
function TESTperformMarkovSteps!(moveidxs,Configs::AbstractMatrix,ConnectedConfs::AbstractArray{T,3},moveWeights::AbstractMatrix,walker_still_running,rng::Random.AbstractRNG) where {T}
    moveidxs = TEST_get_move_idxs(moveidxs,moveWeights, rng)

    # for α in axes(Configs,2)
    AK.foraxes(Configs,2) do α
        if walker_still_running[α]
            move = moveidxs[α]
            Configs[:, α] .= view(ConnectedConfs,:, move, α)
        end
    end
    return nothing
end
function _test_contTime(WE,H,logψ,Hilbert,parallelization)
    return continuos_time_propagation!(WE, H, logψ, Hilbert, 0.1, 0., parallelization)
end

function step_times!(τ_step, τleft, ξ_alpha, el_alpha, Hxx_alpha, logw_alpha,walker_still_running)
    @. τ_step = min(τleft, log(1 - ξ_alpha) ./ (el_alpha - Hxx_alpha))
    @. τleft -= τ_step
    @. walker_still_running = τleft > 0
    @. logw_alpha += -τ_step * el_alpha
end

function step_times_kernel!(τ_step, τleft, ξ_alpha, el_alpha, Hxx_alpha, logw_alpha,walker_still_running)
    @inbounds AK.foreachindex(τ_step) do alpha
        if walker_still_running[alpha]
            τ_step[alpha] = min(τleft[alpha], log(1 - ξ_alpha[alpha]) / (el_alpha[alpha] - Hxx_alpha[alpha]))
            τleft[alpha] -= τ_step[alpha]
            walker_still_running[alpha] = τleft[alpha] > 0
            logw_alpha[alpha] += -τ_step[alpha] * el_alpha[alpha]
        end
    end
end

function continuos_time_propagation_walker!(WE::ArrayWalkerEnsemble, α::Int, H::AbstractSignFreeOperator, logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, dτ::Real, w_avg_estimate::Real,RNG::Random.AbstractRNG)

    Configs = getConfigs(WE)
    moveWeights = getMoveWeights(WE)
    
    AK.foreachindex(Configs,2) do α
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
                error("Infinite propagat
                ion time encountered. Check for too large values in guiding wavefunction or its Buffers!")
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

function _continuos_time_propagation!(WE::ArrayWalkerEnsemble, H::AbstractSignFreeOperator, logψ, Hilbert::AbstractHilbertSpace, dτ::Real, w_avg_estimate::Real, parallelization::KernelParallelization, RNG::Random.AbstractRNG = Random.default_rng())

    MemoryType = getMemoryType(parallelization)

    Configs = getConfigs(WE)
    ConnConfigs = getOffdiag!(WE, H.moves)
    
    moveWeights = getMoveWeights(WE)
    Nwalkers = size(Configs, 2)
    # getMoveWeights!(moveWeights, WE, H, params)
    logw_alpha =  MemoryType{Float64}(undef, Nwalkers) .= 0.0
    Hxx_alpha = MemoryType{Float64}(undef, Nwalkers) .= 0.0
    τ_step =  MemoryType{Float64}(undef, Nwalkers) .= 0.0
    τleft = MemoryType{Float64}(undef, Nwalkers) .= dτ
    ξ_alpha =  MemoryType{Float64}(undef, Nwalkers) .= 0.0
    el_alpha =  MemoryType{Float64}(undef, Nwalkers) .= 0.0
    walker_still_running =  MemoryType{Bool}(undef, Nwalkers) .= true
    moveIdxs = MemoryType{Int}(undef, Nwalkers)
    # walker_still_running_cpu =  ones(Bool, Nwalkers)
    Num_running_walkers = Nwalkers

    Hxx = get_diagonal(H)
    map_function!(Hxx_alpha, Configs, Hxx)

    get_markov_weights!(WE, H, logψ, Hilbert)
    e_local_alpha!(el_alpha, Hxx_alpha, moveWeights)

    i = 0
    a1 = 0.
    a2 = 0.
    a3 = 0.
    a4 = 0.
    # while i < 24
    while AK.any(identity,walker_still_running)
    # while any(walker_still_running)
        i += 1
        Random.rand!(RNG, ξ_alpha)
        a1 += @elapsed step_times_kernel!(τ_step, τleft, ξ_alpha, el_alpha, Hxx_alpha, logw_alpha, walker_still_running)
        # copyto!(walker_still_running_cpu, walker_still_running)
        # if any(isinf,τleft)
        #     # @info "" τ_step el_alpha Hxx_alpha τleft maximum(moveWeights)
        #     error("Infinite propagation time encountered. Check for too large values in guiding wavefunction or its Buffers!")
        # end
        # if !all(walker_still_running)
        #     Num_running_walkers = count(identity, walker_still_running)
        #     @info "Number of still running walkers: $Num_running_walkers"
        # end
        # performMarkovSteps!(Configs, ConnConfigs, moveWeights, walker_still_running, RNG)
        a2 += @elapsed TESTperformMarkovSteps!(moveIdxs,Configs, ConnConfigs, moveWeights, walker_still_running, RNG)
        a3 += @elapsed get_markov_weights!(WE, H, logψ, Hilbert)
        a4 += @elapsed begin 
            map_function!(Hxx_alpha, Configs, Hxx)
            e_local_alpha!(el_alpha, Hxx_alpha, moveWeights)
        end
    end
    print(i,", ", a1*1000,", ", a2*1000,", ", a3*1000,", ", a4*1000,"\n")
    # @info "" i 
    WalkerWeights = getWalkerWeights(WE)
    localEnergies = getLocalEnergies(WE)
    # copyto!(WalkerWeights, exp.(logw_alpha .- dτ * w_avg_estimate))
    # copyto!(localEnergies, el_alpha)
    # @. WalkerWeights = exp(logw_alpha - dτ * w_avg_estimate)
    # @. localEnergies = el_alpha

    return localEnergies
end


function reconfigurateWalkers!(Walkers::ArrayWalkerEnsemble, reconfiguration::MinimalReconfiguration, rng::Random.AbstractRNG)
    reconfigurationList = reconfiguration.reconfigurationList
    reconfiguration_buffer = reconfiguration.reconfigurationBuffer

    Nw = NWalkers(Walkers)
    WalkerWeights = getWalkerWeights(Walkers)
    cumsum!(reconfiguration_buffer, WalkerWeights)
    wTotal = sum(WalkerWeights)
    @. reconfiguration_buffer /= wTotal

    for α in eachindex(Walkers)
        ξα = rand(rng)
        zα = (α + ξα - 1) / Nw
        α′ = searchsortedfirst(reconfiguration_buffer, zα)
        reconfigurationList[α] = α′
    end
    minimizeReconfiguration!(reconfigurationList)

    return nothing
end