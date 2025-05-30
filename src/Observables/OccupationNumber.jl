struct OccupationNumber{T<:Real} <: AbstractObservable
    xBuffer::Vector{T}
end
Base.copy(O::OccupationNumber) = OccupationNumber(copy(O.xBuffer))
obs(O::OccupationNumber) = O.xBuffer
function (O::OccupationNumber)(out,config)
    LoopVectorization.@turbo out .= config
end