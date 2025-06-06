struct OccupationNumber{T<:Real} <: AbstractObservable
    xBuffer::Vector{T}
end
OccupationNumber(Nsites) = OccupationNumber(zeros(Nsites))

Base.copy(O::OccupationNumber) = OccupationNumber(copy(O.xBuffer))
@inline obs(O::OccupationNumber) = O.xBuffer
@inline function (O::OccupationNumber)(out,config::AbstractArray)
    LoopVectorization.@turbo out .= config
end
@inline function (O::OccupationNumber)(out,config::BosonConfig)
    pConf = parent(config)
    O(out, pConf)
end