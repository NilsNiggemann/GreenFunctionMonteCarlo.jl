struct ObservableAccumulator{ObsType<:AbstractObservable,T_high<:AbstractFloat,T_low<:Real} <: AbstractObserver
    BasicAcc::BasicAccumulator{T_high}
    ObsFunc_buffer::Vector{ObsType}
    Obs_Buffers::CircularArrays.CircularArray{T_low, 3, Array{T_low, 3}}
    Obs_numerator::Matrix{T_high}
    Obs_denominator::Vector{T_high}
end

function ObservableAccumulator(filename,Observable::AbstractObservable,m_proj::Integer,NWalkers::Integer,NThreads::Integer)
    p_proj = 2m_proj

    Obs_out = obs(Observable)
    NumObs = length(Obs_out)

    ObsFunc_buffer = [copy(Obs_out) for _ in 1:NThreads]

    Obs_Buffers = CircularArrays.CircularArray(zeros(eltype(Obs_out),NumObs,NWalkers,m_proj))

    Obs_numerator = maybe_MMap_array(filename,"Obs_numerator",Float64,(NumObs,m_proj,))
    Obs_denominator = maybe_MMap_array(filename,"Obs_denominator",Float64,(m_proj,))

    BasicAcc = BasicAccumulator(filename,m_proj,NWalkers)
    ObsAcc = ObservableAccumulator(BasicAcc,ObsFunc_buffer,Obs_Buffers,Obs_numerator,Obs_denominator)
    return ObsAcc
end
ObservableAccumulator(Observable::AbstractObservable,m_proj::Integer,NWalkers::Integer,NThreads::Integer) = ObservableAccumulator(nothing,Observable::AbstractObservable,m_proj::Integer,NWalkers::Integer,NThreads::Integer)

function saveObservables_before!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Hxx = get_diagonal(H)
    energies = Observables.energies
    TotalWeights = Observables.TotalWeights

    update_energies_TotalWeights!(energies,TotalWeights,i,Walkers,Hxx)

    Gnps = Observables.Gnps
    en_numerator = Observables.en_numerator
    en_denominator = Observables.en_denominator

    updateGnp!(Gnps,TotalWeights,i)
    NSites = length(getConfig(Walkers,1))
    getEnergy_step!(en_numerator,en_denominator,Gnps,energies,i,NSites)
    return nothing
end

function saveObservables_after!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Observables.reconfigurationTable[:,i] .= get_reconfigurationList(reconfiguration)
    return nothing
end
