struct MinimalReconfiguration <: AbstractReconfigurationScheme 
    reconfigurationList::Vector{Int}
    reconfigurationBuffer::Vector{Float64}
end
MinimalReconfiguration(Nw::Int) = MinimalReconfiguration(zeros(Int,Nw),zeros(Nw))
get_reconfigurationList(reconf::MinimalReconfiguration) = reconf.reconfigurationList

"""Performs an efficient reconfiguration of walkers. This reconfiguration will not remove walkers if they all have the same weight, which increases the efficiency as more walkers can contribute to the average.

Matteo Calandra Buonaura and Sandro Sorella
Phys. Rev. B 57, 11446 (1998)

# Arguments
- `Walkers::AbstractWalkerEnsemble`: The ensemble of walkers to be reconfigured.
- `reconfigurationList`: A list of indices that will be reconfigured.
- `rng::Random.AbstractRNG`: The random number generator to be used.
"""
function reconfigurateWalkers!(Walkers::AbstractWalkerEnsemble,reconfiguration::MinimalReconfiguration,rng::Random.AbstractRNG)
    reconfigurationList = reconfiguration.reconfigurationList
    reconfiguration_buffer = reconfiguration.reconfigurationBuffer

    Nw = NWalkers(Walkers)
    WalkerWeights = getWalkerWeights(Walkers)
    reconfiguration_buffer = cumsum!(reconfiguration_buffer,WalkerWeights) 
    wTotal = sum(WalkerWeights)
    reconfiguration_buffer ./= wTotal
    for α in eachindex(Walkers)
        ξα = rand(rng)
        zα = (α + ξα - 1)/Nw
        α´ = searchsortedfirst(reconfiguration_buffer,zα)
        reconfigurationList[α] = α´
    end
    minimizeReconfiguration!(reconfigurationList)
    for (α,α´) in enumerate(reconfigurationList)
        if α´ != α
            getConfig(Walkers,α) .= getConfig(Walkers,α´)

            BuffA = getBuffer(Walkers,α)
            BuffB = getBuffer(Walkers,α´)
            setBuffer!(BuffA,BuffB)
        end
    end
end

"""given a list of reconfiguration indices, minimizes the number of reconfigurations by swapping elements in the list. Each walker that survives a reconfiguration step remains unchanged while walkers that are killed get assigned to a new index."""
function minimizeReconfiguration!(list)
    N = length(list)
    index_map = Dict(α′ => α for (α, α′) in enumerate(list))

    for α in 1:N
        α′ = list[α]
        α′ == α && continue

        otherIndex = get(index_map, α, 0)
        iszero(otherIndex) && continue

        list[α], list[otherIndex] = list[otherIndex], list[α]

        index_map[list[otherIndex]] = otherIndex
        index_map[list[α]] = α
    end
    return list
end

function swapIndices!(list,i,j)
    list[i],list[j] = list[j],list[i]
    return list
end