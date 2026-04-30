struct MinimalReconfiguration <: AbstractReconfigurationScheme 
    reconfigurationList::Vector{Int}
    reconfigurationBuffer::Vector{Float64}
    cloned_walkers::Vector{Int}
    dead_walkers::Vector{Int} 
    dead_walker_Set::Set{Int}
end
function MinimalReconfiguration(Nw::Int)
    reconfigurationList = zeros(Int,Nw)
    reconfigurationBuffer = zeros(Nw)
    cloned_walkers = Int[]
    sizehint!(cloned_walkers, Nw)
    dead_walkers = Int[]
    sizehint!(dead_walkers, Nw)
    dead_walker_Set = Set(1:Nw)
    sizehint!(dead_walker_Set, Nw)
    return MinimalReconfiguration(reconfigurationList, reconfigurationBuffer, cloned_walkers, dead_walkers, dead_walker_Set)
end
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
    Polyester.@batch for α in eachindex(Walkers)
        ξα = rand(rng)
        zα = (α + ξα - 1)/Nw
        α´ = searchsortedfirst(reconfiguration_buffer,zα)
        reconfigurationList[α] = α´
    end

    (;cloned_walkers, dead_walkers,dead_walker_Set) = reconfiguration

    minimizeReconfiguration!(reconfigurationList, cloned_walkers, dead_walkers,dead_walker_Set)

    Polyester.@batch for i in eachindex(cloned_walkers,dead_walkers)
        α = dead_walkers[i]
        α´ = cloned_walkers[i]
        getConfig(Walkers,α) .= getConfig(Walkers,α´)

        BuffA = getBuffer(Walkers,α)
        BuffB = getBuffer(Walkers,α´)
        setBuffer!(BuffA,BuffB)
    end
end

"""
    minimizeReconfiguration!(list)


given a list of reconfiguration indices, minimizes the number of reconfigurations by swapping elements in the list. Each walker that survives a reconfiguration step remains unchanged while walkers that are killed get assigned to a new index.

# Arguments
- `list`: A collection (e.g., an array) that will be reconfigured in-place.
"""
function minimizeReconfiguration!(list, cloned_walkers::Vector{Int}, dead_walkers::Vector{Int}, dead_walker_Set::Set{Int})
    N = length(list)

    empty!(dead_walker_Set)
    empty!(dead_walkers)
    for i in 1:N
        push!(dead_walker_Set, i)
    end
    empty!(cloned_walkers)
    
    for α′ in list
        if α′ ∉ dead_walker_Set
            push!(cloned_walkers, α′)
        end
        delete!(dead_walker_Set, α′)
    end
    for α in dead_walker_Set
        push!(dead_walkers, α)
    end
    sort!(cloned_walkers)
    sort!(dead_walkers)
    return cloned_walkers, dead_walkers
end

function swapIndices!(list,i,j)
    list[i],list[j] = list[j],list[i]
    return list
end

struct NoReconfiguration <: AbstractReconfigurationScheme end
reconfigurateWalkers!(Walkers::AbstractWalkerEnsemble,::NoReconfiguration,rng) = nothing