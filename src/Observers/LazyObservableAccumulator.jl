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
    ObsFunc::ObsType
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
    LazyObservableAccumulator(filename,conf<:AbstractConfig , Observable::AbstractObservable, BasicAcc::BasicAccumulator, m_proj::Integer, NWalkers::Integer, NThreads::Integer; Obs_Name = _type_stripped(Observable))

Constructs an `LazyObservableAccumulator` for accumulating measurements of a given observable during a Monte Carlo simulation.

# Arguments
- `filename`: The path to the file where accumulated data will be stored. This argument can be in which case the accumulator will store the result only in memory.
- `conf<:AbstractConfig`: The configuration object used for the simulation.
- `Observable::AbstractObservable`: The observable to be measured and accumulated.
- `BasicAcc::BasicAccumulator`: The basic accumulator object used for storing intermediate results.
- `m_proj::Integer`: The projection quantum number or index relevant to the observable.
- `NWalkers::Integer`: The number of walkers used in the simulation.
- `NThreads::Integer`: The number of threads to be used for parallel accumulation.
- `Obs_Name`: (optional) The name of the observable, defaults to the name of the Observable struct.

# Returns
An `LazyObservableAccumulator` object configured for the specified observable and simulation parameters.

# See also
- [`ObservableAccumulator`](@ref)
"""
function LazyObservableAccumulator(filename,conf::AbstractConfig,Observable::AbstractObservable,BasicAcc::BasicAccumulator,m_values::AbstractVector{Int},NWalkers::Integer,NThreads::Integer;Obs_Name = _type_stripped(Observable))
    Obs_out = obs(Observable)
    NumObs = length(Obs_out)

    num_bins = get_num_bins(BasicAcc)
    ObsFunc = copy(Observable)
    mvals = collect(m_values)


    @assert all(>=(0), mvals) "All m_values must be non-negative integers"
    m_proj = maximum(mvals)+1 # account for 0 projection 
    num_proj = length(mvals)
    @assert m_proj <= projection_order(BasicAcc) "Projection order m_proj=$(m_proj) must be less than the projection order of the BasicAccumulator =$(projection_order(BasicAcc))"
    Obs_Buffers = CircularArrays.CircularArray(zeros(eltype(conf),NWalkers,length(conf),m_proj))
    
    Obs_numerators = maybe_MMap_array(filename,"$(Obs_Name)_numerator",Float64,(NumObs,num_proj,num_bins))
    Obs_denominators = maybe_MMap_array(filename,"$(Obs_Name)_denominator",Float64,(num_proj,num_bins))

    maybe_write_array(filename,"$(Obs_Name)_m_values",mvals)

    ObsAcc = LazyObservableAccumulator(BasicAcc,ObsFunc,Obs_Buffers,Obs_numerators,Obs_denominators,mvals)
    return ObsAcc
end

LazyObservableAccumulator(conf,Observable::AbstractObservable,BasicAcc::BasicAccumulator,m_proj,NWalkers::Integer,NThreads::Integer) = LazyObservableAccumulator(nothing,conf,Observable,BasicAcc,m_proj,NWalkers,NThreads)

function saveObservables_before!(Observables::LazyObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    # saveObservables_before!(Observables.BasicAcc,i,Walkers,H,reconfiguration)
    return nothing
end

function _fill_conf_buffers!(Observables::LazyObservableAccumulator,i,Walkers::AbstractWalkerEnsemble)
    (;Obs_Buffers) = Observables

    Obs_buff_arr = parent(Obs_Buffers)
    i_wrapped = mod1(i,lastindex(Obs_Buffers,3))
    Base.@boundscheck checkbounds(Obs_buff_arr,eachindex(Walkers),eachindex(getConfig(Walkers,1)),i_wrapped)

    Threads.@threads for α in eachindex(Walkers)
        conf = getConfig(Walkers,α)
        @inbounds Obs_buff_arr[α,:,i_wrapped] .= conf
    end

    return
end


function saveObservables_after!(Observables::LazyObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Lazy_Obs_Acc_projection!(Observables,i,Walkers)
    return nothing
end

function Lazy_Obs_Acc_projection!(Observables::LazyObservableAccumulator,n,Walkers::AbstractWalkerEnsemble)
    
    (;Obs_numerators,Obs_denominators,Obs_Buffers) = Observables
    (;PopulationMatrix,Gnps,reconfigurationTable) = Observables.BasicAcc
    
    _fill_conf_buffers!(Observables,n,Walkers)    
    m_max = projection_order(Observables)

    getPopulationMatrix!(PopulationMatrix,reconfigurationTable,n,m_max)
    Nw = length(eachindex(Walkers))
    Obs_Buffers_arr = parent(Obs_Buffers)
    m_values = Observables.m_values
    
    bin_index = get_bin_index(n,Observables.BasicAcc)

    Base.@boundscheck checkbounds(Obs_Buffers_arr,1:Nw,:,:)
    Base.@boundscheck checkbounds(Obs_numerators,:,:,bin_index)
    PopulationMatrix_parent = parent(PopulationMatrix)
    Nw⁻¹ = 1/Nw
    ObsFunc = Observables.ObsFunc
    Threads.@threads for m_index in eachindex(m_values)
        m = m_values[m_index]
        Gnp = Gnps[n,1+2m]
        Gnp == 0 && continue
        Obs_denominators[m_index,bin_index] += Gnp
        n_m_wrapped = mod1(n-m,lastindex(Obs_Buffers,3))
        m_index_wrapped = mod1(m_index,lastindex(PopulationMatrix_parent,2))

        WalkerPopulations = @view PopulationMatrix_parent[:,m_index_wrapped]
        walker_confs = @view Obs_Buffers_arr[:,:,n_m_wrapped]

        Threads.@threads for obs_idx in axes(Obs_numerators,1)

            obs_avg = average_obs_walkers(ObsFunc,obs_idx,walker_confs,WalkerPopulations)

            Obs_numerators[obs_idx,m_index,bin_index] += obs_avg * Nw⁻¹*Gnp
        end
    end
end

"""
    average_obs_walkers(ObsFunc::AbstractObservable, obs_idx::Integer, walker_confs::AbstractMatrix{ConfType}, WalkerPopulations::AbstractVector{<:Integer}) where ConfType
Compute the average value of a observable with index `obs_idx` across a population of walkers, given their configurations and populations.
# Info
This function needs to be implemented for each custom observable type. An efficient implementation is crucial for the performance if the number of observables is large.

# Arguments
- `ObsFunc::AbstractObservable`: The observable function that computes the observable value for a given configuration.
- `obs_idx::Integer`: The index of the observable to be averaged (if the observable returns multiple values, this specifies which one to average).
- `walker_confs::AbstractMatrix{ConfType}`: A matrix Nw × L where Nw is the number of walkers and L is the size of the configuration, containing the configurations of all walkers.
- `WalkerPopulations::AbstractVector{<:Integer}`: A vector containing the population of walkers after a number of reconfigurations, providing the weight the contribution of each configuration to the average.
# Returns
The average value of the specified observable across the walker population, weighted by their populations.

See also
- [`OccupationNumber`](@ref)
- [`SpinCorrelations`](@ref)
"""
average_obs_walkers(ObsFunc,obs_idx,walker_confs,WalkerPopulations) = error("not implemented for custom observable of type $(typeof(ObsFunc))")