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
    Obs_numerator::Matrix{T_high}
    Obs_denominator::Vector{T_high}
end
function reset_accumulator!(Observables::ObservableAccumulator)
    reset_accumulator!(Observables.BasicAcc)
    set_zero!(Observables.Obs_numerator)
    set_zero!(Observables.Obs_denominator)
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
function ObservableAccumulator(filename,Observable::AbstractObservable,BasicAcc::BasicAccumulator,m_proj::Integer,NWalkers::Integer,NThreads::Integer; Obs_Name = _type_stripped(Observable))
    p_proj = 2m_proj
    Obs_out = obs(Observable)
    NumObs = length(Obs_out)

    ObsFunc_buffer = [copy(Observable) for _ in 1:NThreads]
    Obs_Buffers = CircularArrays.CircularArray(zeros(eltype(Obs_out),NumObs,NWalkers,m_proj))

    Obs_numerator = maybe_MMap_array(filename,"$(Obs_Name)_numerator",Float64,(NumObs,m_proj,))
    Obs_denominator = maybe_MMap_array(filename,"$(Obs_Name)_denominator",Float64,(m_proj,))

    ObsAcc = ObservableAccumulator(BasicAcc,ObsFunc_buffer,Obs_Buffers,Obs_numerator,Obs_denominator)
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
    # Obs_Buffers[:,α,i] .= obs_val #todo: allow LoopVectorization for CircularArrays?
    Obs_buff_arr = parent(Obs_Buffers)
    i_wrapped = mod1(i,lastindex(Obs_Buffers,3))
    Base.@boundscheck checkbounds(Obs_buff_arr,1,α,i_wrapped)

    LoopVectorization.@turbo Obs_buff_arr[:,α,i_wrapped] .= obs_val
end

function saveObservables_after!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    # saveObservables_after!(Observables.BasicAcc,i,Walkers,H,reconfiguration)
    Obs_Acc_projection!(Observables,i,Walkers)
    return nothing
end

function Obs_Acc_projection!(Observables::ObservableAccumulator,n,Walkers::AbstractWalkerEnsemble)

    (;Obs_numerator,Obs_denominator,Obs_Buffers) = Observables
    (;PopulationMatrix,Gnps,reconfigurationTable) = Observables.BasicAcc
    
    compute_ObsAccumBuffers!(Observables,n,Walkers)    
    m_max = length(Obs_denominator) - 1

    getPopulationMatrix!(PopulationMatrix,reconfigurationTable,n,m_max)
    Nw = length(eachindex(Walkers))
    Obs_Buffers_arr = parent(Obs_Buffers)

    Nw⁻¹ = 1/Nw
    m_values = 0:m_max
    Threads.@threads for m_index in eachindex(m_values)
        m = m_values[m_index]
        Gnp = Gnps[n,1+2m]
        Gnp == 0 && continue
        # @info "" n m Gnp
        Obs_denominator[m_index] += Gnp
        # Obs_denominator[m_index] += Gnp*Nw
        for α in 1:Nw

            mult = PopulationMatrix[α,m_index]
            mult == 0 && continue
            mult *= Nw⁻¹*Gnp
            n_m_wrapped = mod1(n-m,lastindex(Obs_Buffers,3))
            # O = @view Obs_Buffers_arr[:,α,n_m_wrapped]
            LoopVectorization.@turbo for i in axes(Obs_numerator,1)
                Obs_numerator[i,m_index] += Obs_Buffers_arr[i,α,n_m_wrapped]*mult
            end
        end
    end
end
