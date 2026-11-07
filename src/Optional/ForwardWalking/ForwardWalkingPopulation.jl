# Internal mechanics for one "in-flight" forward-walking sub-population, seeded from the live main
# population by applying an AbstractOffdiagonalObservable once, then continued via the ordinary
# (unmodified) propagateWalkers!/reconfigurateWalkers! primitives for `mProj` further steps. Not
# part of the public API (see ForwardWalkingAccumulator.jl, which owns the estimator math).
#
# Walker-level reweighting: each of the Nsub sub-walkers is seeded with its OWN observable weight
# (`seed_weights[β]`), and `origin[α]` tracks, exactly (via the same discrete-duplication ancestry
# `MinimalReconfiguration` already produces, mirroring `getPopulationMatrix!`'s backward-population
# trick but run forward from the seed), which of the Nsub seed walkers the walker CURRENTLY at slot
# α descends from. This is what lets the estimator (`ForwardWalkingAccumulator`'s `Num`/`Denom`,
# via `origin_weighted_seed_value`) actually correlate a seed configuration's observable value with
# how well its own descendants survive -- collapsing straight to a single population-averaged
# `seed_weight` scalar per seed event (the original design) discards exactly that walker-level
# signal, since population-averaging over Nsub self-averages it away before any reweighting can use
# it (see examples/MWE_forward_walking_bias.jl for the diagnostic that found this).
mutable struct ForwardWalkingSlot{WE<:AbstractWalkerEnsemble}
    Walkers::WE
    reconfiguration::MinimalReconfiguration
    active::Bool
    n_seed::Int      # main-loop step index this slot was seeded at
    step::Int        # 0 (just seeded) .. mProj
    seed_weights::Vector{Float64}  # seed_weights[β]: observable weight of seed walker β
    origin::Vector{Int}            # origin[α]: which seed walker (1..Nsub) current walker α descends from
    origin_scratch::Vector{Int}    # reusable buffer for _propagate_origin!
    cum_weight::Float64            # population-level growth factor (shared across all Nsub walkers)
end

function allocate_forward_walking_slot(conf, logψ::AbstractGuidingFunction, Nsub::Integer, H::AbstractSignFreeOperator)
    Walkers = allocate_walkerEnsemble(conf, logψ, Nsub, H)
    reconfiguration = MinimalReconfiguration(Nsub)
    seed_weights = zeros(Nsub)
    origin = collect(1:Nsub)
    origin_scratch = collect(1:Nsub)
    return ForwardWalkingSlot(Walkers, reconfiguration, false, 0, 0, seed_weights, origin, origin_scratch, 1.0)
end

# resolve_moves=true support: one "seed group" is n_moves independent ForwardWalkingSlots (one per
# move), seeded/advanced/freed in lockstep -- see ForwardWalkingAccumulator.jl's type docstring.
allocate_forward_walking_group(conf, logψ::AbstractGuidingFunction, Nsub::Integer, H::AbstractSignFreeOperator, n_moves::Integer) =
    [allocate_forward_walking_slot(conf, logψ, Nsub, H) for _ in 1:n_moves]

# Seeds `slot` from the CURRENTLY LIVE (pre-reconfiguration) main population `Walkers`, applying the
# off-diagonal observable once per sub-walker. Direct analog of SpiderWebModel's
# `initialize_forward_walking!`, but reading from the live population instead of a stored snapshot.
function seed_from_main_population!(slot::ForwardWalkingSlot, Walkers::AbstractWalkerEnsemble, Observable::AbstractOffdiagonalObservable,
                                     logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, n::Integer,
                                     weights_buf::AbstractVector, rng::Random.AbstractRNG)
    Nsub = NWalkers(slot.Walkers)
    @assert Nsub == NWalkers(Walkers) "ForwardWalkingAccumulator requires Nsub == NWalkers(main population) for v1 (1:1 seed, no resampling); got Nsub=$Nsub, NWalkers(main)=$(NWalkers(Walkers))"

    for β in 1:Nsub
        Config = getConfig(slot.Walkers, β)
        Config .= getConfig(Walkers, β)
        Buffer = getBuffer(slot.Walkers, β)
        compute_GWF_buffer!(Buffer, logψ, Config)
        get_observable_weights!(weights_buf, Config, Observable, logψ, Hilbert, Buffer)
        move, w = sample_and_apply_observable!(Config, weights_buf, Observable, rng)
        move === nothing || post_move_affect!(Buffer, Config, move, logψ)
        slot.seed_weights[β] = w
    end
    _reweight_seed!(slot.seed_weights, getWalkerWeights(Walkers))

    slot.origin .= eachindex(slot.origin)
    slot.n_seed = n
    slot.step = 0
    slot.cum_weight = 1.0
    slot.active = true
    return slot
end

# resolve_moves=true variant: seeds all n_moves slots in `group` at once, one per move, by applying
# EVERY move (not stochastically sampling one) to its own copy of each main-population walker. Unlike
# `sample_and_apply_observable!`'s single-sample sum estimator, this gives each slot the exact
# per-move mixed-estimator seed weight `ψ(x'_i)/ψ(x)` for every sub-walker, since nothing here is
# randomly chosen.
function seed_from_main_population!(group::AbstractVector{<:ForwardWalkingSlot}, Walkers::AbstractWalkerEnsemble, Observable::AbstractOffdiagonalObservable,
                                     logψ::AbstractGuidingFunction, Hilbert::AbstractHilbertSpace, n::Integer, rng::Random.AbstractRNG)
    n_mv = n_moves(Observable)
    @assert length(group) == n_mv "resolve_moves seed group size ($(length(group))) must match n_moves(Observable) ($n_mv)"

    Nsub = NWalkers(group[1].Walkers)
    @assert Nsub == NWalkers(Walkers) "ForwardWalkingAccumulator requires Nsub == NWalkers(main population) for v1 (1:1 seed, no resampling); got Nsub=$Nsub, NWalkers(main)=$(NWalkers(Walkers))"

    for β in 1:Nsub
        x = getConfig(Walkers, β)
        for i in 1:n_mv
            slot = group[i]
            Config = getConfig(slot.Walkers, β)
            Config .= x
            Buffer = getBuffer(slot.Walkers, β)
            compute_GWF_buffer!(Buffer, logψ, Config)
            pre_move_affect!(Buffer, Config, logψ)
            move = get_move(Observable, i)
            Δ = log_psi_diff(Config, move, logψ, Buffer, Hilbert)
            ψratio = exp(Δ)
            w = observable_weight(Observable, i, ψratio, Config)
            apply!(Config, move)
            post_move_affect!(Buffer, Config, move, logψ)
            slot.seed_weights[β] = w
        end
    end
    foreach(slot -> _reweight_seed!(slot.seed_weights, getWalkerWeights(Walkers)), group)

    for slot in group
        slot.origin .= eachindex(slot.origin)
        slot.n_seed = n
        slot.step = 0
        slot.cum_weight = 1.0
        slot.active = true
    end
    return group
end

# The seed is taken after propagation but before reconfiguration, so walker β still carries its
# step weight w_β (the weight the reconfiguration would resample by, and that the diagonal mixed
# estimator `getLocalEnergyWalkers_before` uses). Neither `Num` (cum_weight starts after the seed
# step) nor `Denom` (Gnps window starts at n_seed+1) contains it, so it must enter the seed value:
# omitting it biases ⟨O⟩ by ≈ -cov(w, O_loc)/w̄ (∝ dτ, zero for an exact guide).
_reweight_seed!(seed_weights, walker_weights) = (seed_weights .*= walker_weights ./ Statistics.mean(walker_weights))

# Propagates ancestry labels through one reconfiguration step: the walker now at slot α is a copy
# of whatever was at pre-reconfiguration slot `reconfigurationList[α]`, so it inherits THAT slot's
# origin label. Mirrors the same ancestry-chain idea `getPopulationMatrix!`
# (src/Observables/computeObservables.jl) uses for the main population's diagonal observables, run
# forward from the seed instead of backward from "now".
function _propagate_origin!(origin::Vector{Int}, scratch::Vector{Int}, reconfigurationList::Vector{Int})
    scratch .= origin
    @inbounds for α in eachindex(origin)
        origin[α] = scratch[reconfigurationList[α]]
    end
    return origin
end

# Advances `slot` by exactly one propagation+reconfiguration step, reusing the ordinary propagator
# and reconfiguration scheme completely unmodified. Updates `slot.cum_weight` (the shared
# population-level growth factor) and `slot.origin` (exact per-walker seed ancestry).
function advance_slot!(slot::ForwardWalkingSlot, propagator::AbstractPropagator, logψ::AbstractGuidingFunction,
                        H::AbstractSignFreeOperator, Hilbert::AbstractHilbertSpace, w_avg_estimate::Real,
                        parallelization::AbstractParallelizationScheme, rng::Random.AbstractRNG)
    # propagateWalkers!/reconfigurateWalkers! take a Vector of per-task RNGs (to support
    # MultiThreaded/BatchMultiThreaded), matching how runGFMC!'s own main loop derives RNGs from a
    # single rng via duplicate_rng -- mirror that here rather than passing FW's scalar rng directly.
    RNGs = duplicate_rng(rng, num_tasks(parallelization))
    propagateWalkers!(slot.Walkers, H, logψ, Hilbert, propagator, parallelization, RNGs)
    branch_weight = Statistics.mean(getWalkerWeights(slot.Walkers))
    slot.cum_weight *= branch_weight / w_avg_estimate
    reconfigurateWalkers!(slot.Walkers, slot.reconfiguration, parallelization, RNGs)
    _propagate_origin!(slot.origin, slot.origin_scratch, get_reconfigurationList(slot.reconfiguration))
    slot.step += 1
    return slot
end

# The forward-walking estimator's per-slot contribution to Num at the slot's current projection
# order: the EXACT descendant-weighted average of the seed observable values (see the type
# docstring). Σ_α seed_weights[origin[α]] groups by seed walker exactly like
# Σ_β count[β]*seed_weights[β] would (count[β] = number of current walkers descended from seed
# walker β, which sums to Nsub), without needing to build that histogram explicitly.
function origin_weighted_seed_value(slot::ForwardWalkingSlot)
    (; seed_weights, origin) = slot
    acc = 0.0
    @inbounds for α in eachindex(origin)
        acc += seed_weights[origin[α]]
    end
    return acc / length(origin)
end

function _find_free_slot(slots::AbstractVector{<:ForwardWalkingSlot})
    for slot in slots
        slot.active || return slot
    end
    error("No free forward-walking slot available; this indicates mProj/cadence sizing is inconsistent with the slot pool size.")
end

# resolve_moves=true variant: every slot in a group is seeded/advanced/freed in lockstep, so checking
# the first slot's `active` flag suffices.
function _find_free_slot(slot_groups::AbstractVector{<:AbstractVector{<:ForwardWalkingSlot}})
    for group in slot_groups
        group[1].active || return group
    end
    error("No free forward-walking slot group available; this indicates mProj/cadence sizing is inconsistent with the slot pool size.")
end
