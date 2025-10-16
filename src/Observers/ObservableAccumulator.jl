"""
    ObservableAccumulator{ObsType<:AbstractObservable, T_high<:AbstractFloat, T_low<:Real}

An accumulator for observables in Monte Carlo simulations. This struct is designed to collect and process measurements of a given observable type during the simulation, supporting both high-precision (`T_high`) and lower-precision (`T_low`) data types.

# Type Parameters
- `ObsType<:AbstractObservable`: The type of observable being accumulated.
- `T_high<:AbstractFloat`: The floating-point type used for high-precision accumulation (e.g., `Float64`).
- `T_low<:Real`: The type used for lower-precision to speed-up calculations (typically `Float32`).

# Description
`ObservableAccumulator` is typically used as part of the observer pattern in Monte Carlo simulations, where it collects measurements of observables at each step and provides methods for statistical analysis, such as computing averages and variances.


# See Also
- [`BasicAccumulator`](@ref)
"""
struct ObservableAccumulator{ObsType<:AbstractObservable,T_high<:AbstractFloat,T_low<:Real} <: AbstractObserver
    BasicAcc::BasicAccumulator{T_high}
    ObsFunc_buffer::Vector{ObsType}
    Obs_Buffers::CircularArrays.CircularArray{T_low, 3, Array{T_low, 3}}
    Obs_numerators::Array{T_high,3}
    Obs_denominators::Matrix{T_high}
end
projection_order(Observables::ObservableAccumulator) = size(Observables.Obs_denominators,1) - 1

function reset_accumulator!(Observables::ObservableAccumulator)
    set_zero!(Observables.Obs_numerators)
    set_zero!(Observables.Obs_denominators)
    return Observables
end
_type_stripped(::T) where T = nameof(T)
"""
    ObservableAccumulator(filename, Observable::AbstractObservable, BasicAcc::BasicAccumulator, m_proj::Integer, NWalkers::Integer, NThreads::Integer; Obs_Name = _type_stripped(Observable))

Constructs an `ObservableAccumulator` for accumulating measurements of a given observable during a Monte Carlo simulation.

# Arguments
- `filename`: The path to the file where accumulated data will be stored. This argument can be in which case the accumulator will store the result only in memory.
- `Observable::AbstractObservable`: The observable to be measured and accumulated.
- `BasicAcc::BasicAccumulator`: The basic accumulator object used for storing intermediate results.
- `m_proj::Integer`: The projection quantum number or index relevant to the observable.
- `NWalkers::Integer`: The number of walkers used in the simulation.
- `NThreads::Integer`: The number of threads to be used for parallel accumulation.
- `Obs_Name`: (optional) The name of the observable, defaults to the name of the Observable struct.

# Returns
An `ObservableAccumulator` object configured for the specified observable and simulation parameters.
"""
function ObservableAccumulator(filename,Observable::AbstractObservable,BasicAcc::BasicAccumulator,m_proj::Integer,NWalkers::Integer,NThreads::Integer;Obs_Name = _type_stripped(Observable))
    Obs_out = obs(Observable)
    NumObs = length(Obs_out)

    num_bins = get_num_bins(BasicAcc)

    ObsFunc_buffer = [copy(Observable) for _ in 1:NThreads]
    Obs_Buffers = CircularArrays.CircularArray(zeros(eltype(Obs_out),NumObs,NWalkers,m_proj))

    Obs_numerators = maybe_MMap_array(filename,"$(Obs_Name)_numerator",Float64,(NumObs,m_proj,num_bins))
    Obs_denominators = maybe_MMap_array(filename,"$(Obs_Name)_denominator",Float64,(m_proj,num_bins))

    ObsAcc = ObservableAccumulator(BasicAcc,ObsFunc_buffer,Obs_Buffers,Obs_numerators,Obs_denominators)
    return ObsAcc
end

ObservableAccumulator(Observable::AbstractObservable,BasicAcc::BasicAccumulator,m_proj::Integer,NWalkers::Integer,NThreads::Integer) = ObservableAccumulator(nothing,Observable::AbstractObservable,BasicAcc,m_proj::Integer,NWalkers::Integer,NThreads::Integer)

function saveObservables_before!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    # saveObservables_before!(Observables.BasicAcc,i,Walkers,H,reconfiguration)
    return nothing
end

function compute_ObsAccumBuffers!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble)
    numThreads = length(Observables.ObsFunc_buffer)
    if numThreads == 1
        _compute_ObsAccumBuffers_singlethreaded!(Observables,i,Walkers)
    else
        _compute_ObsAccumBuffers_multithreaded!(Observables,i,Walkers)
    end
    return
end

function _compute_ObsAccumBuffers_singlethreaded!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble)
    (;ObsFunc_buffer,Obs_Buffers) = Observables

    ObsFunc! = ObsFunc_buffer[begin]

    for α in eachindex(Walkers)
        conf = getConfig(Walkers,α)
        _kernel_compute_ObsAccumBuffers!(Obs_Buffers,conf,α,i,ObsFunc!)
    end
    return
end

function _compute_ObsAccumBuffers_multithreaded!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble)
    (;ObsFunc_buffer,Obs_Buffers) = Observables

    batches = ChunkSplitters.chunks(eachindex(Walkers), n = length(ObsFunc_buffer),split = ChunkSplitters.RoundRobin())

    @sync for (i_chunk, αinds) in enumerate(batches)
        ObsFunc! = ObsFunc_buffer[i_chunk]

        Threads.@spawn for α in αinds
            conf = getConfig(Walkers,α)
            _kernel_compute_ObsAccumBuffers!(Obs_Buffers,conf,α,i,ObsFunc!)
        end
    end
    return
end

Base.@propagate_inbounds function _kernel_compute_ObsAccumBuffers!(Obs_Buffers,conf,α,i,ObsFunc!)
    obs_val = obs(ObsFunc!)
    ObsFunc!(obs_val,conf)
    Obs_buff_arr = parent(Obs_Buffers)
    i_wrapped = mod1(i,lastindex(Obs_Buffers,3))
    Base.@boundscheck checkbounds(Obs_buff_arr,:,α,i_wrapped)

    LoopVectorization.@turbo Obs_buff_arr[:,α,i_wrapped] .= obs_val
end

function saveObservables_after!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    # saveObservables_after!(Observables.BasicAcc,i,Walkers,H,reconfiguration)
    Obs_Acc_projection!(Observables,i,Walkers)
    return nothing
end

function Obs_Acc_projection!(Observables::ObservableAccumulator,n,Walkers::AbstractWalkerEnsemble)

    (;Obs_numerators,Obs_denominators,Obs_Buffers) = Observables
    (;PopulationMatrix,Gnps,reconfigurationTable) = Observables.BasicAcc
    
    compute_ObsAccumBuffers!(Observables,n,Walkers)    
    m_max = projection_order(Observables)

    getPopulationMatrix!(PopulationMatrix,reconfigurationTable,n,m_max)
    Nw = length(eachindex(Walkers))
    Obs_Buffers_arr = parent(Obs_Buffers)

    Nw⁻¹ = 1/Nw
    m_values = 0:m_max
    
    bin_index = get_bin_index(n,Observables.BasicAcc)

    Base.@boundscheck checkbounds(Obs_Buffers_arr,axes(Obs_numerators,1),1:Nw,:)
    Base.@boundscheck checkbounds(m_values,axes(Obs_numerators,2))
    Base.@boundscheck checkbounds(Obs_numerators,:,:,bin_index)
    PopulationMatrix_parent = parent(PopulationMatrix)
    Nw⁻¹ = 1/Nw
    Threads.@threads for m_index in axes(Obs_numerators,2)
        m = m_values[m_index]
        Gnp = Gnps[n,1+2m]
        Gnp == 0 && continue
        # @info "" n m Gnp
        Obs_denominators[m_index] += Gnp
        # Obs_denominator[m_index] += Gnp*Nw
        for α in 1:Nw

            mult = PopulationMatrix[α,m_index]
            mult == 0 && continue
            mult *= Nw⁻¹*Gnp
            n_m_wrapped = mod1(n-m,lastindex(Obs_Buffers,3))
            # O = @view Obs_Buffers_arr[:,α,n_m_wrapped]
            LoopVectorization.@turbo for i in axes(Obs_numerators,1)
                Obs_numerators[i,m_index] += Obs_Buffers_arr[i,α,n_m_wrapped]*mult
            end
        end
    end
end

@views function get_obs_from_accumulator(Obs::Union{ObservableAccumulator,NamedTuple},bin_indices::AbstractVector)
    # Obs_num_slices = Obs.Obs_numerators[:,:,bin_indices]
    # Obs_denom_slices = Obs.Obs_denominators[:,bin_indices]

    Normalization = Statistics.mean(Obs.Obs_denominators[1,:])

    Obs_num = zeros(eltype(Obs.Obs_numerators),size(Obs.Obs_numerators,1),size(Obs.Obs_numerators,2))
    Obs_denom = zeros(eltype(Obs.Obs_denominators),size(Obs.Obs_denominators,1))

    for bin_idx in bin_indices
        Obs_num .+= Obs.Obs_numerators[:,:,bin_idx] ./ Normalization
        Obs_denom .+= Obs.Obs_denominators[:,bin_idx] ./ Normalization 
    end
    Obs_num ./= Obs_denom'
    return Obs_num
end
get_obs_from_accumulator(Observables::Union{ObservableAccumulator,NamedTuple}) = [get_obs_from_accumulator(Observables,idx:idx) for idx in axes(Observables.Obs_denominators,2)]
function get_obs_from_accumulator_bunching(Observables::Union{ObservableAccumulator,NamedTuple},n_bunch::Integer;kwargs...)
    chunks = ChunkSplitters.chunks(axes(Observables.Obs_denominators,2), size = n_bunch, split = ChunkSplitters.Consecutive();kwargs...)
    return [
        get_obs_from_accumulator(Observables,chunk)
        for chunk in chunks
    ]
end