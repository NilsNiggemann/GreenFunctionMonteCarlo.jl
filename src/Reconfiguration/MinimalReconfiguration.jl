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
- `RNGs::Vector{<:Random.AbstractRNG}`: A vector of random number generators to be used.
"""
function reconfigurateWalkers!(Walkers::AbstractWalkerEnsemble,reconfiguration::MinimalReconfiguration,parallelization::AbstractParallelizationScheme,RNGs::Vector{<:Random.AbstractRNG})
    reconfigurationList = reconfiguration.reconfigurationList
    reconfiguration_buffer = reconfiguration.reconfigurationBuffer
    Nw⁻¹ = 1. /NWalkers(Walkers)
    WalkerWeights = getWalkerWeights(Walkers)
    reconfiguration_buffer = cumsum!(reconfiguration_buffer,WalkerWeights) 
    wTotal = sum(WalkerWeights)
    reconfiguration_buffer ./= wTotal

    _make_reconfigurationList!(reconfigurationList,reconfiguration_buffer, Nw⁻¹, RNGs, parallelization)

    (;cloned_walkers, dead_walkers,dead_walker_Set) = reconfiguration

    minimizeReconfiguration!(reconfigurationList, cloned_walkers, dead_walkers,dead_walker_Set)

    _replace_walkers!(Walkers,cloned_walkers,dead_walkers,parallelization)
end
function _replace_walker!(Walkers,α,α´)
    getConfig(Walkers,α) .= getConfig(Walkers,α´)

    BuffA = getBuffer(Walkers,α)
    BuffB = getBuffer(Walkers,α´)
    setBuffer!(BuffA,BuffB)
end
function _replace_walkers!(Walkers,cloned_walkers,dead_walkers,::AbstractParallelizationScheme)
    for i in eachindex(cloned_walkers,dead_walkers)
        α = dead_walkers[i]
        α´ = cloned_walkers[i]
        _replace_walker!(Walkers,α,α´)
    end
end
function _replace_walkers!(Walkers,cloned_walkers,dead_walkers,::BatchMultiThreaded)
    Polyester.@batch per=thread for i in eachindex(cloned_walkers,dead_walkers)
        α = dead_walkers[i]
        α´ = cloned_walkers[i]
        _replace_walker!(Walkers,α,α´)
    end
end

function _make_reconfigurationList!(reconfigurationList,reconfiguration_buffer, Nw⁻¹, RNGs, ::AbstractParallelizationScheme)
    for α in eachindex(reconfigurationList)
        reconfigurationList[α] = _sample_reconf(reconfiguration_buffer,α, Nw⁻¹, first(RNGs))
    end
end

function _make_reconfigurationList!(reconfigurationList,reconfiguration_buffer, Nw⁻¹, RNGs, parallelization::BatchMultiThreaded)
    Nw = length(reconfigurationList)
    minbatch = clamp(Nw ÷ num_tasks(parallelization), 1, Nw)

    Polyester.@batch per=thread minbatch=minbatch for α in eachindex(reconfigurationList)
        task_idx = polyester_get_task_local_idx(α, minbatch)
        rng = RNGs[task_idx]
        reconfigurationList[α] = _sample_reconf(reconfiguration_buffer,α, Nw⁻¹, rng)
    end
end

function _sample_reconf(reconfiguration_buffer,α, Nw⁻¹, rng::Random.AbstractRNG)
    ξα = rand(rng)
    zα = (α + ξα - 1)*Nw⁻¹
    α´ = searchsortedfirst(reconfiguration_buffer,zα)
    return α´
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

    # `list` (== reconfigurationList, later saved verbatim into reconfigurationTable) must keep
    # recording each walker's TRUE post-reconfiguration source. `_replace_walkers!` only physically
    # overwrites `dead_walkers` slots (from the independently-sorted `cloned_walkers` pairing below,
    # which no longer matches the original per-walker sample in `list`), leaving every other slot's
    # own pre-reconfiguration value untouched -- so `list` must be rewritten to match exactly that:
    # self-reference for untouched slots, the actual clone source for overwritten ones. Backward
    # population tracing (getPopulationMatrix!) depends on `list` being an accurate historical
    # record, not just on the walker ensemble ending up with the statistically-correct population.
    for α in 1:N
        list[α] = α
    end
    for i in eachindex(dead_walkers, cloned_walkers)
        list[dead_walkers[i]] = cloned_walkers[i]
    end

    return cloned_walkers, dead_walkers
end

function swapIndices!(list,i,j)
    list[i],list[j] = list[j],list[i]
    return list
end

struct NoReconfiguration <: AbstractReconfigurationScheme end
reconfigurateWalkers!(Walkers::AbstractWalkerEnsemble,::NoReconfiguration,parallelization::AbstractParallelizationScheme,rng) = nothing