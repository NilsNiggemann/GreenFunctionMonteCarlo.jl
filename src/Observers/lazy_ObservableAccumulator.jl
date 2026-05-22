"""
    LazyObservableAccumulator{ObsType<:AbstractObservable, T_high<:AbstractFloat, T_low<:Real}

An accumulator for observables in Monte Carlo simulations. This struct is designed to collect and process measurements of a given observable type during the simulation.

# Type Parameters
- `ObsType<:AbstractObservable`: The type of observable being accumulated.
- `T_high<:AbstractFloat`: The floating-point type used for high-precision accumulation (e.g., `Float64`).
- `T_low<:Real`: The type used for buffering configs (typically `Int8`).

# Description
`LazyObservableAccumulator` is typically used as part of the observer pattern in Monte Carlo simulations, where it collects measurements of observables at each step and provides methods for statistical analysis, such as computing averages and variances.


# See Also
- [`BasicAccumulator`](@ref)
- [`ObservableAccumulator`](@ref)
"""
struct LazyObservableAccumulator{ObsType<:AbstractObservable,T_high<:AbstractFloat,T_low<:Real} <: AbstractObserver
    BasicAcc::BasicAccumulator{T_high}
    ObsFunc_buffer::Vector{ObsType}
    Obs_Buffers::CircularArrays.CircularArray{T_low, 3, Array{T_low, 3}}
    Obs_numerators::Array{T_high,3}
    Obs_denominators::Matrix{T_high}
    m_values::Vector{Int}
end
projection_order(Observables::LazyObservableAccumulator) = size(Observables.Obs_denominators,1) - 1

function reset_accumulator!(Observables::LazyObservableAccumulator)
    set_zero!(Observables.Obs_numerators)
    set_zero!(Observables.Obs_denominators)
    return Observables
end
"""
    LazyObservableAccumulator(filename, Observable::AbstractObservable, BasicAcc::BasicAccumulator, m_proj::Integer, NWalkers::Integer, NThreads::Integer; Obs_Name = _type_stripped(Observable))

Constructs an `LazyObservableAccumulator` for accumulating measurements of a given observable during a Monte Carlo simulation.

# Arguments
- `filename`: The path to the file where accumulated data will be stored. This argument can be in which case the accumulator will store the result only in memory.
- `Observable::AbstractObservable`: The observable to be measured and accumulated.
- `BasicAcc::BasicAccumulator`: The basic accumulator object used for storing intermediate results.
- `m_proj::Integer`: The projection quantum number or index relevant to the observable.
- `NWalkers::Integer`: The number of walkers used in the simulation.
- `NThreads::Integer`: The number of threads to be used for parallel accumulation.
- `Obs_Name`: (optional) The name of the observable, defaults to the name of the Observable struct.

# Returns
An `LazyObservableAccumulator` object configured for the specified observable and simulation parameters.
"""
function LazyObservableAccumulator(filename,conf,Observable::AbstractObservable,BasicAcc::BasicAccumulator,m_values::AbstractVector{Int},NWalkers::Integer,NThreads::Integer;Obs_Name = _type_stripped(Observable))
    Obs_out = obs(Observable)
    NumObs = length(Obs_out)

    num_bins = get_num_bins(BasicAcc)
    ObsFunc_buffer = [copy(Observable) for _ in 1:NThreads]
    mvals = collect(m_values)


    @assert all(>=(0), mvals) "All m_values must be non-negative integers"
    m_proj = maximum(mvals)+1 # account for 0 projection 
    num_proj = length(mvals)
    @assert m_proj <= projection_order(BasicAcc) "Projection order m_proj=$(m_proj) must be less than the projection order of the BasicAccumulator =$(projection_order(BasicAcc))"
    Obs_Buffers = CircularArrays.CircularArray(zeros(eltype(conf),length(conf),NWalkers,m_proj))
    
    Obs_numerators = maybe_MMap_array(filename,"$(Obs_Name)_numerator",Float64,(NumObs,num_proj,num_bins))
    Obs_denominators = maybe_MMap_array(filename,"$(Obs_Name)_denominator",Float64,(num_proj,num_bins))

    maybe_write_array(filename,"$(Obs_Name)_m_values",mvals)

    ObsAcc = LazyObservableAccumulator(BasicAcc,ObsFunc_buffer,Obs_Buffers,Obs_numerators,Obs_denominators,mvals)
    return ObsAcc
end

LazyObservableAccumulator(Observable::AbstractObservable,BasicAcc::BasicAccumulator,m_proj,NWalkers::Integer,NThreads::Integer) = LazyObservableAccumulator(nothing,Observable::AbstractObservable,BasicAcc,m_proj::Integer,NWalkers::Integer,NThreads::Integer)

function saveObservables_before!(Observables::LazyObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    # saveObservables_before!(Observables.BasicAcc,i,Walkers,H,reconfiguration)
    return nothing
end

function compute_ObsAccumBuffers!(Observables::LazyObservableAccumulator,i,Walkers::AbstractWalkerEnsemble)
    numThreads = length(Observables.ObsFunc_buffer)
    if numThreads == 1
        _compute_ObsAccumBuffers_singlethreaded!(Observables,i,Walkers)
    else
        _compute_ObsAccumBuffers_multithreaded!(Observables,i,Walkers)
    end
    return
end

function _compute_ObsAccumBuffers_singlethreaded!(Observables::LazyObservableAccumulator,i,Walkers::AbstractWalkerEnsemble)
    (;ObsFunc_buffer,Obs_Buffers) = Observables

    ObsFunc! = ObsFunc_buffer[begin]

    for α in eachindex(Walkers)
        conf = getConfig(Walkers,α)
        _kernel_compute_ObsAccumBuffers!(Obs_Buffers,conf,α,i,ObsFunc!)
    end
    return
end

function _compute_ObsAccumBuffers_multithreaded!(Observables::LazyObservableAccumulator,i,Walkers::AbstractWalkerEnsemble)
    (;ObsFunc_buffer,Obs_Buffers) = Observables

    batches = ChunkSplitters.chunks(eachindex(Walkers), n = length(ObsFunc_buffer),split = ChunkSplitters.RoundRobin())

    # @sync for (i_chunk, αinds) in enumerate(batches)
    #     Threads.@spawn for α in αinds
    #         conf = getConfig(Walkers,α)
    #         _kernel_compute_ObsAccumBuffers!(Obs_Buffers,conf,α,i)
    #     end
    # end

    Polyester.@batch for α in eachindex(Walkers)
        conf = getConfig(Walkers,α)
        _kernel_compute_ObsAccumBuffers!(Obs_Buffers,conf,α,i)
    end
    return
end

@inline Base.@propagate_inbounds function _kernel_compute_ObsAccumBuffers!(Obs_Buffers,conf,α,i)
    Obs_buff_arr = parent(Obs_Buffers)
    i_wrapped = mod1(i,lastindex(Obs_Buffers,3))
    Base.@boundscheck checkbounds(Obs_buff_arr,:,α,i_wrapped)

    Obs_buff_arr[:,α,i_wrapped] .= conf
    # LoopVectorization.@turbo Obs_buff_arr[:,α,i_wrapped] .= conf
end

function saveObservables_after!(Observables::LazyObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    # saveObservables_after!(Observables.BasicAcc,i,Walkers,H,reconfiguration)
    Obs_Acc_projection!(Observables,i,Walkers)
    return nothing
end

function Obs_Acc_projection!(Observables::LazyObservableAccumulator,n,Walkers::AbstractWalkerEnsemble)

    (;Obs_numerators,Obs_denominators,Obs_Buffers) = Observables
    (;PopulationMatrix,Gnps,reconfigurationTable) = Observables.BasicAcc
    
    nThreads = length(Observables.ObsFunc_buffer)
    compute_ObsAccumBuffers!(Observables,n,Walkers)    
    m_max = projection_order(Observables)

    getPopulationMatrix!(PopulationMatrix,reconfigurationTable,n,m_max)
    Nw = length(eachindex(Walkers))
    Obs_Buffers_arr = parent(Obs_Buffers)
    m_values = Observables.m_values
    
    bin_index = get_bin_index(n,Observables.BasicAcc)

    Base.@boundscheck checkbounds(Obs_Buffers_arr,:,1:Nw,:)
    # Base.@boundscheck checkbounds(axes(Obs_numerators,2),m_values)
    Base.@boundscheck checkbounds(Obs_numerators,:,:,bin_index)
    PopulationMatrix_parent = parent(PopulationMatrix)
    Nw⁻¹ = 1/Nw
    # batches = ChunkSplitters.chunks(eachindex(m_values), n = nThreads,split = ChunkSplitters.RoundRobin())
    # error(collect.(batches))
    i_chunk = 1        
    ObsFunc! = Observables.ObsFunc_buffer[i_chunk]
    obs_val = obs(ObsFunc!)

    Threads.@threads for m_index in eachindex(m_values)
        m = m_values[m_index]
        Gnp = Gnps[n,1+2m]
        Gnp == 0 && continue
        Obs_denominators[m_index,bin_index] += Gnp
        n_m_wrapped = mod1(n-m,lastindex(Obs_Buffers,3))
        m_index_wrapped = mod1(m_index,lastindex(PopulationMatrix_parent,2))

        Threads.@threads for i in axes(Obs_numerators,1)

            Obs_num_i_m_b = zero(eltype(Obs_numerators))
            
            for α in 1:Nw
                mult = PopulationMatrix_parent[α,m_index_wrapped]
                mult == 0 && continue
                mult *= Nw⁻¹*Gnp
                x_am = @view Obs_Buffers_arr[:,α,n_m_wrapped]
                Obs_num_i_m_b += ObsFunc!(x_am,i) *mult
            end
            Obs_numerators[i,m_index,bin_index] += Obs_num_i_m_b
        end
    end
end

# function Obs_Acc_projection!(Observables::LazyObservableAccumulator,n,Walkers::AbstractWalkerEnsemble)

#     (;Obs_numerators,Obs_denominators,Obs_Buffers) = Observables
#     (;PopulationMatrix,Gnps,reconfigurationTable) = Observables.BasicAcc
    
#     nThreads = length(Observables.ObsFunc_buffer)
#     compute_ObsAccumBuffers!(Observables,n,Walkers)    
#     m_max = projection_order(Observables)

#     getPopulationMatrix!(PopulationMatrix,reconfigurationTable,n,m_max)
#     Nw = length(eachindex(Walkers))
#     Obs_Buffers_arr = parent(Obs_Buffers)
#     m_values = Observables.m_values
    
#     bin_index = get_bin_index(n,Observables.BasicAcc)

#     Base.@boundscheck checkbounds(Obs_Buffers_arr,:,1:Nw,:)
#     # Base.@boundscheck checkbounds(axes(Obs_numerators,2),m_values)
#     Base.@boundscheck checkbounds(Obs_numerators,:,:,bin_index)
#     PopulationMatrix_parent = parent(PopulationMatrix)
#     Nw⁻¹ = 1/Nw
#     batches = ChunkSplitters.chunks(eachindex(m_values), n = nThreads,split = ChunkSplitters.RoundRobin())
#     # error(collect.(batches))
#     @sync for (i_chunk, m_batch) in enumerate(batches)
        
#         Threads.@spawn begin
#             ObsFunc! = Observables.ObsFunc_buffer[i_chunk]
#             obs_val = obs(ObsFunc!)
#             for m_index in m_batch
#                 m = m_values[m_index]
#                 Gnp = Gnps[n,1+2m]
#                 Gnp == 0 && continue
#                 Obs_denominators[m_index,bin_index] += Gnp
#                 n_m_wrapped = mod1(n-m,lastindex(Obs_Buffers,3))
#                 m_index_wrapped = mod1(m_index,lastindex(PopulationMatrix_parent,2))
#                 for α in 1:Nw
#                     mult = PopulationMatrix_parent[α,m_index_wrapped]
#                     mult == 0 && continue
#                     mult *= Nw⁻¹*Gnp
#                     @views ObsFunc!(obs_val, Obs_Buffers_arr[:,α,n_m_wrapped])

#                     LoopVectorization.@tturbo for i in axes(Obs_numerators,1)
#                         Obs_numerators[i,m_index,bin_index] += obs_val[i]*mult
#                     end
#                 end
#             end
#         end
#     end
# end

# @views function get_obs_from_accumulator(Obs::Union{LazyObservableAccumulator,NamedTuple},bin_indices::AbstractVector)
#     # Obs_num_slices = Obs.Obs_numerators[:,:,bin_indices]
#     # Obs_denom_slices = Obs.Obs_denominators[:,bin_indices]

#     Normalization = Statistics.mean(Obs.Obs_denominators[1,:])

#     Obs_num = zeros(eltype(Obs.Obs_numerators),size(Obs.Obs_numerators,1),size(Obs.Obs_numerators,2))
#     Obs_denom = zeros(eltype(Obs.Obs_denominators),size(Obs.Obs_denominators,1))

#     for bin_idx in bin_indices
#         Obs_num .+= Obs.Obs_numerators[:,:,bin_idx] ./ Normalization
#         Obs_denom .+= Obs.Obs_denominators[:,bin_idx] ./ Normalization 
#     end
#     Obs_num ./= Obs_denom'
#     return Obs_num
# end
# get_obs_from_accumulator(Observables::Union{LazyObservableAccumulator,NamedTuple}) = [get_obs_from_accumulator(Observables,idx:idx) for idx in axes(Observables.Obs_denominators,2)]
# function get_obs_from_accumulator_bunching(Observables::Union{LazyObservableAccumulator,NamedTuple},n_bunch::Integer;kwargs...)
#     chunks = ChunkSplitters.chunks(axes(Observables.Obs_denominators,2), size = n_bunch, split = ChunkSplitters.Consecutive();kwargs...)
#     return [
#         get_obs_from_accumulator(Observables,chunk)
#         for chunk in chunks
#     ]
# end