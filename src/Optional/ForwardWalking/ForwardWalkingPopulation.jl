# Internal mechanics for one "in-flight" forward-walking sub-population, seeded from the live main
# population by applying an AbstractOffdiagonalObservable once, then continued via the ordinary
# (unmodified) propagateWalkers!/reconfigurateWalkers! primitives for `mProj` further steps. Not
# part of the public API (see ForwardWalkingAccumulator.jl, which owns the estimator math).
mutable struct ForwardWalkingSlot{WE<:AbstractWalkerEnsemble}
    Walkers::WE
    reconfiguration::MinimalReconfiguration
    active::Bool
    n_seed::Int      # main-loop step index this slot was seeded at
    step::Int        # 0 (just seeded) .. mProj
    seed_weight::Float64
    cum_weight::Float64
end

function allocate_forward_walking_slot(conf, logψ::AbstractGuidingFunction, Nsub::Integer, H::AbstractSignFreeOperator)
    Walkers = allocate_walkerEnsemble(conf, logψ, Nsub, H)
    reconfiguration = MinimalReconfiguration(Nsub)
    return ForwardWalkingSlot(Walkers, reconfiguration, false, 0, 0, 0.0, 1.0)
end

# Seeds `slot` from the CURRENTLY LIVE (pre-reconfiguration) main population `Walkers`, applying the
# off-diagonal observable once per sub-walker. Direct analog of SpiderWebModel's
# `initialize_forward_walking!`, but reading from the live population instead of a stored snapshot.
function seed_from_main_population!(slot::ForwardWalkingSlot, Walkers::AbstractWalkerEnsemble, Observable::AbstractOffdiagonalObservable,
                                     logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, n::Integer,
                                     weights_buf::AbstractVector, rng::Random.AbstractRNG)
    Nsub = NWalkers(slot.Walkers)
    @assert Nsub == NWalkers(Walkers) "ForwardWalkingAccumulator requires Nsub == NWalkers(main population) for v1 (1:1 seed, no resampling); got Nsub=$Nsub, NWalkers(main)=$(NWalkers(Walkers))"

    total_weight = 0.0
    for β in 1:Nsub
        Config = getConfig(slot.Walkers, β)
        Config .= getConfig(Walkers, β)
        Buffer = getBuffer(slot.Walkers, β)
        compute_GWF_buffer!(Buffer, logψ, Config)
        get_observable_weights!(weights_buf, Config, Observable, logψ, Hilbert, Buffer)
        move, w = sample_and_apply_observable!(Config, weights_buf, Observable, rng)
        move === nothing || post_move_affect!(Buffer, Config, move, logψ)
        total_weight += w
    end

    slot.seed_weight = total_weight / Nsub
    slot.n_seed = n
    slot.step = 0
    slot.cum_weight = 1.0
    slot.active = true
    return slot
end

# Advances `slot` by exactly one propagation+reconfiguration step, reusing the ordinary propagator
# and reconfiguration scheme completely unmodified. Returns the sub-population's mean branching
# weight for this step (folded into `slot.cum_weight`).
function advance_slot!(slot::ForwardWalkingSlot, propagator::AbstractPropagator, logψ::AbstractGuidingFunction,
                        H::AbstractSignFreeOperator, Hilbert::AbstractHilbertSpace, w_avg_estimate::Real,
                        parallelization::AbstractParallelizationScheme, rng::Random.AbstractRNG)
    propagateWalkers!(slot.Walkers, H, logψ, Hilbert, propagator, parallelization, rng)
    branch_weight = Statistics.mean(getWalkerWeights(slot.Walkers))
    slot.cum_weight *= branch_weight / w_avg_estimate
    reconfigurateWalkers!(slot.Walkers, slot.reconfiguration, rng)
    slot.step += 1
    return slot
end

function _find_free_slot(slots::AbstractVector{<:ForwardWalkingSlot})
    for slot in slots
        slot.active || return slot
    end
    error("No free forward-walking slot available; this indicates mProj/cadence sizing is inconsistent with the slot pool size.")
end
