struct OccupationNumber{T<:Real} <: AbstractObservable
    xBuffer::Vector{T}
end
OccupationNumber(Nsites) = OccupationNumber(zeros(Float32,Nsites))

Base.copy(O::OccupationNumber) = OccupationNumber(copy(O.xBuffer))
@inline obs(O::OccupationNumber) = O.xBuffer
@inline function (O::OccupationNumber)(out,config::AbstractArray)
    LoopVectorization.@turbo out .= config
end
@inline function (O::OccupationNumber)(out,config::BosonConfig)
    pConf = parent(config)
    O(out, pConf)
end

function average_obs_walkers(O::OccupationNumber,obs_idx::Integer,walker_confs::AbstractMatrix, WalkerPopulations::AbstractVector{<:Integer})
    Obs_num_i_m_b = zero(eltype(obs(O)))

    i = obs_idx
    LoopVectorization.@turbo for α in axes(walker_confs,1)
        mult = WalkerPopulations[α]
        x_i = walker_confs[α,i]

        Obs_num_i_m_b += x_i * mult
    end
    return Obs_num_i_m_b
end
