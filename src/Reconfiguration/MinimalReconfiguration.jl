struct MinimalReconfiguration <: AbstractReconfigurationScheme 
    reconfigurationList::Vector{Int}
    reconfigurationBuffer::Vector{Float64}
    inverseBuffer::Vector{Int}
end
MinimalReconfiguration(Nw::Int) = MinimalReconfiguration(zeros(Int,Nw),zeros(Nw),zeros(Int,Nw))
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
    minimizeReconfiguration!(reconfigurationList, reconfiguration.inverseBuffer)
    for (α,α´) in enumerate(reconfigurationList)
        if α´ != α
            getConfig(Walkers,α) .= getConfig(Walkers,α´)

            BuffA = getBuffer(Walkers,α)
            BuffB = getBuffer(Walkers,α´)
            setBuffer!(BuffA,BuffB)
        end
    end
end

"""
    minimizeReconfiguration!(list, inverse)


given a list of reconfiguration indices, minimizes the number of reconfigurations by swapping elements in the list. Each walker that survives a reconfiguration step remains unchanged while walkers that are killed get assigned to a new index.

`inverse` is a pre-allocated scratch buffer of the same length as `list`.  It
is overwritten during the call and its contents after the call are undefined.
Passing a pre-allocated buffer avoids the heap allocation that would otherwise
occur on every call.

# Arguments
- `list`: A collection (e.g., an array) that will be reconfigured in-place.
- `inverse`: A pre-allocated integer scratch buffer of the same length as `list`.
"""
function minimizeReconfiguration!(list, inverse::AbstractVector{<:Integer})
    N = length(list)
    # Clear the buffer so untargeted positions remain 0
    fill!(inverse, 0)
    # Build inverse permutation: inverse[α′] = α  ⟺  list[α] = α′
    @inbounds for α in 1:N
        inverse[list[α]] = α
    end
    @inbounds for α in 1:N
        α′ = list[α]
        α′ == α && continue

        otherIndex = inverse[α]
        iszero(otherIndex) && continue

        list[α], list[otherIndex] = list[otherIndex], list[α]

        inverse[list[otherIndex]] = otherIndex
        inverse[list[α]] = α
    end
    return list
end

"""
    minimizeReconfiguration!(list)

Convenience single-argument overload that allocates a temporary inverse buffer.
Prefer the two-argument form `minimizeReconfiguration!(list, inverse)` in
performance-critical code to avoid heap allocation.
"""
function minimizeReconfiguration!(list)
    return minimizeReconfiguration!(list, zeros(Int, length(list)))
end

function swapIndices!(list,i,j)
    list[i],list[j] = list[j],list[i]
    return list
end

struct NoReconfiguration <: AbstractReconfigurationScheme end
reconfigurateWalkers!(Walkers::AbstractWalkerEnsemble,::NoReconfiguration,rng) = nothing