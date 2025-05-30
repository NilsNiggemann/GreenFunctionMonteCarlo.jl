test() = println("Runni.")

struct ObservableAccumulator{ObsType<:AbstractObservable,T_high<:AbstractFloat,T_low<:Real} <: AbstractObserver
    BasicAcc::BasicAccumulator{T_high}
    ObsFunc_buffer::Vector{ObsType}
    Obs_Buffers::CircularArrays.CircularArray{T_low, 3, Array{T_low, 3}}
    Obs_numerator::Matrix{T_high}
    Obs_denominator::Vector{T_high}
end

function ObservableAccumulator(filename,Observable::AbstractObservable,BasicAcc::BasicAccumulator,m_proj::Integer,NWalkers::Integer,NThreads::Integer)
    p_proj = 2m_proj
    Obs_out = obs(Observable)
    NumObs = length(Obs_out)

    ObsFunc_buffer = [copy(Observable) for _ in 1:NThreads]
    Obs_Buffers = CircularArrays.CircularArray(zeros(eltype(Obs_out),NumObs,NWalkers,m_proj))

    Obs_numerator = maybe_MMap_array(filename,"Obs_numerator",Float64,(NumObs,m_proj,))
    Obs_denominator = maybe_MMap_array(filename,"Obs_denominator",Float64,(m_proj,))

    ObsAcc = ObservableAccumulator(BasicAcc,ObsFunc_buffer,Obs_Buffers,Obs_numerator,Obs_denominator)
    return ObsAcc
end

ObservableAccumulator(Observable::AbstractObservable,m_proj::Integer,NWalkers::Integer,NThreads::Integer) = ObservableAccumulator(nothing,Observable::AbstractObservable,m_proj::Integer,NWalkers::Integer,NThreads::Integer)

function saveObservables_before!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    # saveObservables_before!(Observables.BasicAcc,i,Walkers,H,reconfiguration)
    return nothing
end

function compute_ObsAccumBuffers!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble)
    (;ObsFunc_buffer,Obs_Buffers) = Observables

    batches = ChunkSplitters.chunks(eachindex(Walkers), n = Threads.nthreads(),split= ChunkSplitters.RoundRobin())

    @sync for (i_chunk, αinds) in enumerate(batches)
        ObsFunc! = ObsFunc_buffer[i_chunk]

        Threads.@spawn for α in αinds
            conf = getConfig(Walkers,α)
            Obs_view = @view Obs_Buffers[:,α,i]

            obs_val = obs(ObsFunc!)
            ObsFunc!(obs_val,conf)
            LoopVectorization.@turbo Obs_view .= obs_val
        end
    end
    return
end

function saveObservables_after!(Observables::ObservableAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    # saveObservables_after!(Observables.BasicAcc,i,Walkers,H,reconfiguration)
    Obs_Acc_projection!(Observables,i,Walkers)
    return nothing
end

function Obs_Acc_projection!(Observables::ObservableAccumulator,n,Walkers::AbstractWalkerEnsemble)

    (;Obs_numerator,Obs_denominator,Obs_Buffers) = Observables
    (;PopulationMatrix,Gnps,reconfigurationTable) = Observables.BasicAcc
    (;PopulationMatrix,Gnps,reconfigurationTable) = Observables.BasicAcc
    
    compute_ObsAccumBuffers!(Observables,n,Walkers)    
    m_max = length(Obs_denominator) - 1

    getPopulationMatrix!(PopulationMatrix,reconfigurationTable,n,m_max)
    Nw = length(eachindex(Walkers))
    
    Nw⁻¹ = 1/Nw

    m_values = 0:m_max
    Threads.@sync for m_index in eachindex(m_values)
        Threads.@spawn begin
            m = m_values[m_index]
            Gnp = Gnps[n,1+2m]
            Obs_denominator[m_index] += Gnp
            # Obs_denominator[m_index] += Gnp*Nw
            @views for α in 1:Nw
                mult = PopulationMatrix[α,m_index]
                mult == 0 && continue
                mult *= Nw⁻¹
                O = Obs_Buffers[:,α,n-m]
                @. Obs_numerator[:,m_index] += O*Gnp*mult
            end
        end
    end
end
