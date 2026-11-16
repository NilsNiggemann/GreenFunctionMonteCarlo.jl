"""
    ForwardWalkingAccumulator{ObsType<:AbstractOffdiagonalObservable} <: AbstractObserver

A streaming accumulator for the equal-time expectation value `⟨O⟩` of an off-diagonal operator `O`
(an [`AbstractOffdiagonalObservable`](@ref)), via forward walking.

Unlike diagonal observables (see [`ObservableAccumulator`](@ref)), an off-diagonal operator changes
the sampled configuration when applied, so its expectation value cannot be measured "for free" and
reweighted using only backward population-ancestry bookkeeping. Instead, periodically (every
`cadence` main-loop steps), a small disposable sub-population is seeded from the *live*, currently
running main population, the operator is applied once to seed it, and it is then evolved forward for
`mProj` further steps using the ordinary propagator and reconfiguration scheme — reusing
`propagateWalkers!`/`reconfigurateWalkers!` completely unmodified. At most `fld(mProj,cadence)+1`
such sub-populations are ever in flight, giving bounded memory with no need to store the full main
chain's configuration history.

At forward-projection depth `p=0` (the seed weight alone, no continuation) this is the standard
biased mixed estimator `⟨ψ_T|O|ψ_0⟩`. As `p` grows toward `mProj`, continuing real dynamics
re-projects the seeded state toward the true ground state — the same bias-removal mechanism as the
existing diagonal `Gnp[n,p]` energy estimator. Concretely, `Num` is the forked sub-population's own
weight growth (`slot.cum_weight`) times its ancestry-weighted seed observable value, while `Denom` is
the *reference* (main, unperturbed) population's own weight growth over the same window, read from
`BasicAcc.Gnps` — reusing the paper's insight (Sec. XI) that the already-running main population is
itself a valid, unbiased continuation to normalize against, so no separate denominator-only chain
needs to be run. `BasicAcc.Gnps[n,p]` is a *backward*-looking product ending at step `n`; since a
forward-walking slot's `step` advances in lockstep with the main loop from the moment it's seeded,
reading it at the *current* step `n` (not the seed step `n_seed`) gives exactly the main population's
growth over `[n_seed+1, n]` — the same window the forked sub-population has itself been advancing
over. See `straight_forward_walking.md` (repo root) for the full derivation.

Crucially, the bias-correction is *walker-level*, not population-level: each of the `Nsub` seed
walkers keeps its own observable weight, and each sub-walker's exact seed ancestry is tracked
forward (via `MinimalReconfiguration`'s discrete-duplication list, the same ancestry mechanism
[`ObservableAccumulator`](@ref)'s `PopulationMatrix` uses backward for diagonal observables) rather
than collapsing the whole sub-population into one population-averaged seed value up front —
population-averaging before reweighting discards exactly the per-walker signal the correction needs
(see `origin_weighted_seed_value` in `ForwardWalkingPopulation.jl`, and
`examples/MWE_forward_walking_bias.jl`, which is what found and confirms the fix for this). Check
convergence of [`get_observable_from_accumulator`](@ref)'s output vs. projection order `p` to
validate that `mProj` is large enough.

# Requirement
`FW.BasicAcc` must be the *same* [`BasicAccumulator`](@ref) instance also included, separately, as
its own observer in the enclosing [`CombinedObserver`](@ref), and it must run *before*
`ForwardWalkingAccumulator` in that `CombinedObserver` — exactly as already required for
[`ObservableAccumulator`](@ref)/[`LazyObservableAccumulator`](@ref). This is now load-bearing for
correctness, not just bin bookkeeping: every call reads `BasicAcc.Gnps[n, ·]` for the *current* step
`n`, which is only valid once `BasicAcc`'s own `saveObservables_before!` has already updated it for
that step.

# Multi-move observables and `resolve_moves`
When `Observable` covers more than one move (see [`offdiagonalObservable`](@ref)), the `resolve_moves`
constructor keyword controls whether the result is their combined *sum* (default) or a genuine
per-move breakdown — see the constructor docstring below.
"""
struct ForwardWalkingAccumulator{ObsType<:AbstractOffdiagonalObservable,SlotType,Prop<:AbstractPropagator,GF<:AbstractGuidingFunction,HS<:AbstractHilbertSpace,PS<:AbstractParallelizationScheme,T_high<:AbstractFloat,NA<:AbstractArray{Float64},RNG<:Random.AbstractRNG,FS<:ForwardWalkingSlot} <: AbstractObserver
    BasicAcc::BasicAccumulator{T_high}
    Observable::ObsType
    slots::Vector{SlotType}
    cadence::Int
    mProj::Int
    propagator::Prop
    logψ::GF
    Hilbert::HS
    parallelization::PS
    weights_bufs::Vector{Vector{Float64}}  # per-thread scratch for seeding
    RNGs::Vector{RNG}                       # per-thread RNGs, allocated once
    pending::Vector{FS}                     # reusable list of slots to advance in the current step
    Num::NA
    Denom::Matrix{Float64}
end

"""
    ForwardWalkingAccumulator(filename, Observable::AbstractOffdiagonalObservable, BasicAcc::BasicAccumulator,
                               conf, logψ::AbstractGuidingFunction, H::AbstractSignFreeOperator, Hilbert::AbstractHilbertSpace,
                               propagator::AbstractPropagator, mProj::Integer, cadence::Integer, Nsub::Integer;
                               parallelization=SingleThreaded(), rng=Random.default_rng(), num_bins=get_num_bins(BasicAcc))

Construct a [`ForwardWalkingAccumulator`](@ref) for measuring `⟨Observable⟩` via streaming forward walking.

# Arguments
- `filename`: HDF5 file to persist `Num`/`Denom` to (mirroring `BasicAccumulator`/`ObservableAccumulator`), or `nothing` for an in-memory-only accumulator.
- `Observable::AbstractOffdiagonalObservable`: the off-diagonal operator to measure.
- `BasicAcc::BasicAccumulator`: supplies both the shared bin bookkeeping (`get_bin_index`/`get_num_bins`) and, via `BasicAcc.Gnps`, the reference-population weight growth used to build `Denom` (see the type docstring). Must also be included, separately, in the enclosing `CombinedObserver`, and must run *before* this accumulator there.
- `conf`: an exemplary configuration, used to allocate each forward-walking sub-population.
- `logψ`, `H`, `Hilbert`, `propagator`: the same guiding function, Hamiltonian, Hilbert space, and propagator used by the main simulation — sub-populations are continued with *ordinary* dynamics after being seeded.
- `mProj::Integer`: number of forward-walking continuation steps (the projection-order sweep goes from `p=0` to `p=mProj`). Must satisfy `mProj < 2*projection_order(BasicAcc)`, so that `BasicAcc.Gnps`'s circular buffer has enough valid columns to cover every continuation depth (asserted at construction).
- `cadence::Integer`: seed a new forward-walking sub-population every `cadence` main-loop steps.
- `Nsub::Integer`: size of each forward-walking sub-population. Must equal the main population's `NWalkers` (v1: 1:1 seeding, no resampling).

# Keyword Arguments
- `parallelization`: parallelization scheme for seeding and advancing sub-populations. Default `SingleThreaded()`. Any other scheme runs seeding over walkers, propagation over all (sub-population, walker) pairs, and reconfiguration over sub-populations as Polyester batches. When several `ForwardWalkingAccumulator`s sit in one `CombinedObserver`, their sub-populations are advanced together in one joint pass (see `advance_forward_walking!`).
- `rng`: random number generator from which the per-thread RNGs used for seeding/advancing sub-populations are derived once, at construction. Note that `saveObservables_before!`/`_after!` do not receive the main loop's RNG, so this accumulator owns its own streams. With a multithreaded `parallelization` the result is not bitwise reproducible (walkers are assigned to threads dynamically), as for the main loop.
- `num_bins`: number of bins for error-bar bunching (default: `BasicAcc`'s own).
- `resolve_moves::Bool`: if `false` (default), `Observable`'s moves are combined into a single
  scalar estimate per projection order — for a single-move `Observable` this is `⟨move⟩` itself; for
  a multi-move `Observable`, each seed event stochastically samples *one* move (categorically, ∝ its
  weight) and the accumulated result is an unbiased estimator of the **sum** of the observable over
  all moves (e.g. a structure factor `Σᵢ⟨Oᵢ⟩`), never a per-move breakdown. If `true`, each seed event
  instead branches **deterministically** into one forward-walking sub-population per move (applying
  every move separately, not sampling one), and `get_observable_from_accumulator`/
  `get_observable_from_accumulator_bunching` return one estimate *per move* (an `(n_moves, mProj+1)`
  matrix per bin) — the exact per-move analogue of running `n_moves` separate single-move
  `ForwardWalkingAccumulator`s, bundled into one object/array. Costs `n_moves` times the sub-population
  memory/compute of the default mode (there is no cheaper way to recover individual move values, since
  the default mode's stochastic sampling only ever estimates their sum).
"""
function ForwardWalkingAccumulator(filename, Observable::AbstractOffdiagonalObservable, BasicAcc::BasicAccumulator,
                                    conf, logψ::AbstractGuidingFunction, H::AbstractSignFreeOperator, Hilbert::AbstractHilbertSpace,
                                    propagator::AbstractPropagator, mProj::Integer, cadence::Integer, Nsub::Integer;
                                    parallelization::AbstractParallelizationScheme = SingleThreaded(),
                                    rng::Random.AbstractRNG = Random.default_rng(),
                                    num_bins::Integer = get_num_bins(BasicAcc),
                                    resolve_moves::Bool = false)

    @assert mProj < 2 * projection_order(BasicAcc) "ForwardWalkingAccumulator requires mProj < 2*projection_order(BasicAcc), got mProj=$mProj, 2*projection_order(BasicAcc)=$(2*projection_order(BasicAcc)); otherwise BasicAcc.Gnps's circular buffer doesn't have enough valid columns to cover every forward-walking continuation depth, silently corrupting Denom."

    n_slots = fld(mProj, cadence) + 1
    n_mv = n_moves(Observable)
    n_scratch = _fw_num_scratch(parallelization)
    weights_bufs = [zeros(n_mv) for _ in 1:n_scratch]
    RNGs = duplicate_rng(rng, n_scratch)

    Denom = maybe_MMap_array(filename, "$(_type_stripped(Observable))_ForwardWalking_denominator", Float64, (mProj + 1, num_bins))

    if resolve_moves
        slots = [allocate_forward_walking_group(conf, logψ, Nsub, H, n_mv) for _ in 1:n_slots]
        Num = maybe_MMap_array(filename, "$(_type_stripped(Observable))_ForwardWalking_numerator", Float64, (n_mv, mProj + 1, num_bins))
    else
        slots = [allocate_forward_walking_slot(conf, logψ, Nsub, H) for _ in 1:n_slots]
        Num = maybe_MMap_array(filename, "$(_type_stripped(Observable))_ForwardWalking_numerator", Float64, (mProj + 1, num_bins))
    end

    pending = similar(_flat_slots(slots), 0)
    return ForwardWalkingAccumulator(BasicAcc, Observable, slots, Int(cadence), Int(mProj), propagator, logψ, Hilbert,
                                      parallelization, weights_bufs, RNGs, pending, Num, Denom)
end

function ForwardWalkingAccumulator(Observable::AbstractOffdiagonalObservable, BasicAcc::BasicAccumulator, args...; kwargs...)
    return ForwardWalkingAccumulator(nothing, Observable, BasicAcc, args...; kwargs...)
end

_flat_slots(slots::AbstractVector{<:ForwardWalkingSlot}) = slots
_flat_slots(groups::AbstractVector{<:AbstractVector{<:ForwardWalkingSlot}}) = reduce(vcat, groups)

# Appends to `out` every active slot of `FW` that has not yet been advanced in this step.
function _collect_pending_slots!(out::AbstractVector, FW::ForwardWalkingAccumulator)
    for slot_or_group in FW.slots
        for slot in _as_group(slot_or_group)
            slot.active && !slot.advanced && push!(out, slot)
        end
    end
    return out
end
_as_group(slot::ForwardWalkingSlot) = (slot,)
_as_group(group::AbstractVector{<:ForwardWalkingSlot}) = group

# FW accumulators that can be advanced in one joint pass with `ref`: same propagation dynamics,
# reference population and slot type.
_can_advance_jointly(ref::ForwardWalkingAccumulator, FW::ForwardWalkingAccumulator) =
    FW.BasicAcc === ref.BasicAcc && FW.propagator == ref.propagator && FW.logψ === ref.logψ && FW.Hilbert === ref.Hilbert && eltype(FW.pending) === eltype(ref.pending)

_forward_walking_accumulators() = ()
_forward_walking_accumulators(FW::ForwardWalkingAccumulator, rest...) = (FW, _forward_walking_accumulators(rest...)...)
_forward_walking_accumulators(C::CombinedObserver, rest...) = (_forward_walking_accumulators(C.Observers...)..., _forward_walking_accumulators(rest...)...)
_forward_walking_accumulators(_, rest...) = _forward_walking_accumulators(rest...)

"""
    advance_forward_walking!(observers, n, H)

Advance, by one step, every active forward-walking sub-population of all
[`ForwardWalkingAccumulator`](@ref)s among `observers` (a tuple; nested `CombinedObserver`s are
searched too) in one joint pass: a single parallel loop over all (sub-population, walker) pairs, then
a single parallel loop over sub-populations for reconfiguration. Called by `CombinedObserver` before
its observers record, so that ~`n_accumulators * (fld(mProj, cadence) + 1)` small parallel regions per
step become two large ones. Accumulators whose propagator, guide, Hilbert space or `BasicAcc` differ
from the first one's are left alone and advance themselves in their own `saveObservables_before!`.
Uses the first accumulator's `parallelization` and RNGs.
"""
advance_forward_walking!(observers::Tuple, n, H::AbstractSignFreeOperator) = _advance_jointly!(_forward_walking_accumulators(observers...), H)
advance_forward_walking!(observers, n, H::AbstractSignFreeOperator) = advance_forward_walking!(Tuple(observers), n, H)

_advance_jointly!(::Tuple{}, H::AbstractSignFreeOperator) = nothing
function _advance_jointly!(FWs::Tuple, H::AbstractSignFreeOperator)
    ref = first(FWs)
    pending = empty!(ref.pending)
    for FW in FWs
        _can_advance_jointly(ref, FW) && _collect_pending_slots!(pending, FW)
    end
    advance_slots!(pending, ref.propagator, ref.logψ, H, ref.Hilbert, ref.BasicAcc.weight_normalization[], ref.parallelization, ref.RNGs)
    foreach(slot -> slot.advanced = true, pending)
    empty!(pending)
    return nothing
end

function saveObservables_before!(FW::ForwardWalkingAccumulator, n, Walkers::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, reconfiguration::AbstractReconfigurationScheme)
    # Slots not already advanced by a joint pass (standalone use, or an incompatible accumulator).
    _advance_jointly!((FW,), H)
    _record_slots!(FW, n)
    if mod(n, FW.cadence) == 0
        slot = _find_free_slot(FW.slots)
        _seed!(slot, FW, Walkers, n)
        n > 1 && _record_seed!(FW, slot, get_bin_index(n, FW.BasicAcc))
    end
    return nothing
end

_seed!(slot::ForwardWalkingSlot, FW::ForwardWalkingAccumulator, Walkers::AbstractWalkerEnsemble, n) =
    seed_from_main_population!(slot, Walkers, FW.Observable, FW.logψ, FW.Hilbert, n, FW.weights_bufs, FW.RNGs, FW.parallelization)
_seed!(group::AbstractVector{<:ForwardWalkingSlot}, FW::ForwardWalkingAccumulator, Walkers::AbstractWalkerEnsemble, n) =
    seed_from_main_population!(group, Walkers, FW.Observable, FW.logψ, FW.Hilbert, n, FW.parallelization)

# Adds the contributions of all active slots (just advanced, at their current depth `step`) to
# Num/Denom, frees slots that reached mProj, and clears the `advanced` flags for the next step.
function _record_slots!(FW::ForwardWalkingAccumulator, n)
    bin_index = get_bin_index(n, FW.BasicAcc)
    for slot_or_group in FW.slots
        group = _as_group(slot_or_group)
        first(group).active || continue
        # No warm-up guard needed here (unlike BasicAccumulator's own getEnergy_step!, which guards
        # `n > p`): a slot's `step` advances by exactly 1 per main-loop call from the moment it's
        # seeded, so the current step is always `n = n_seed + slot.step = n_seed + c - 1`, and since
        # `n_seed >= 1` this gives `n >= c`, which is always within Gnps[n,c]'s freshly-written,
        # valid range (see the type docstring for why reading at the *current* n, not n_seed, is what
        # makes this the correct forward-looking reference-population weight growth).
        c = first(group).step + 1
        _add_num!(FW.Num, slot_or_group, c, bin_index)
        # Denom is the reference-population growth, independent of the observable/move, so it is
        # shared across every move of a resolve_moves group.
        FW.Denom[c, bin_index] += FW.BasicAcc.Gnps[n, c]
        for slot in group
            slot.advanced = false
            slot.step == FW.mProj && (slot.active = false)
        end
    end
    return nothing
end

_add_num!(Num::AbstractMatrix, slot::ForwardWalkingSlot, c, bin_index) = (Num[c, bin_index] += origin_weighted_seed_value(slot) * slot.cum_weight)
function _add_num!(Num::AbstractArray{<:Real,3}, group::AbstractVector{<:ForwardWalkingSlot}, c, bin_index)
    for (i, slot) in enumerate(group)
        Num[i, c, bin_index] += origin_weighted_seed_value(slot) * slot.cum_weight
    end
    return Num
end

function _record_seed!(FW::ForwardWalkingAccumulator, slot::ForwardWalkingSlot, bin_index)
    FW.Num[1, bin_index] += Statistics.mean(slot.seed_weights)
    FW.Denom[1, bin_index] += slot.cum_weight
    return nothing
end
function _record_seed!(FW::ForwardWalkingAccumulator, group::AbstractVector{<:ForwardWalkingSlot}, bin_index)
    for (i, slot) in enumerate(group)
        FW.Num[i, 1, bin_index] += Statistics.mean(slot.seed_weights)
    end
    FW.Denom[1, bin_index] += first(group).cum_weight
    return nothing
end

saveObservables_after!(FW::ForwardWalkingAccumulator, i, Walkers::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, reconfiguration::AbstractReconfigurationScheme) = nothing

"""
    get_observable_from_accumulator(FW::ForwardWalkingAccumulator, bin_indices::AbstractVector)
    get_observable_from_accumulator(FW::ForwardWalkingAccumulator)

Compute `⟨Observable⟩` at every forward-projection order `p = 0, 1, ..., mProj` from the accumulated
`Num`/`Denom` state, optionally restricted to a subset of bins (for error-bar analysis). With no
`bin_indices` argument, returns one estimate per bin.
"""
@views function get_observable_from_accumulator(FW::ForwardWalkingAccumulator, bin_indices::AbstractVector)
    Normalization = Statistics.mean(FW.Denom[1, :])

    num = zeros(eltype(FW.Num), size(FW.Num, 1))
    denom = zeros(eltype(FW.Denom), size(FW.Denom, 1))

    for bin_idx in bin_indices
        num .+= FW.Num[:, bin_idx] ./ Normalization
        denom .+= FW.Denom[:, bin_idx] ./ Normalization
    end
    num ./= denom
    return num
end
get_observable_from_accumulator(FW::ForwardWalkingAccumulator) = [get_observable_from_accumulator(FW, idx:idx) for idx in axes(FW.Denom, 2)]

# resolve_moves=true variant: FW.Num is (n_moves, mProj+1, num_bins), FW.Denom stays (mProj+1,
# num_bins) shared across moves (see the type docstring) -- exactly the shape `get_obs_from_accumulator`
# (ObservableAccumulator.jl) already handles generically, so just forward to it. More specific than
# the `bin_indices::AbstractVector` method above, so it takes priority for resolve_moves accumulators;
# the zero-arg `get_observable_from_accumulator(FW)` and `get_observable_from_accumulator_bunching`
# above already delegate to this one via dispatch, so no further overloads are needed.
get_observable_from_accumulator(FW::ForwardWalkingAccumulator{ObsType,SlotType} where {ObsType<:AbstractOffdiagonalObservable,SlotType<:AbstractVector{<:ForwardWalkingSlot}}, bin_indices::AbstractVector) = get_obs_from_accumulator(FW.Num, FW.Denom, bin_indices)

"""
    get_observable_from_accumulator_bunching(FW::ForwardWalkingAccumulator, n_bunch::Integer; kwargs...)

Like [`get_observable_from_accumulator`](@ref), but groups the accumulated bins into `n_bunch` groups
first (for statistical error analysis), mirroring [`get_energy_from_accumulator_bunching`](@ref).
"""
function get_observable_from_accumulator_bunching(FW::ForwardWalkingAccumulator, n_bunch::Integer; kwargs...)
    chunks = ChunkSplitters.chunks(axes(FW.Denom, 2), size = n_bunch, split = ChunkSplitters.Consecutive(); kwargs...)
    return [get_observable_from_accumulator(FW, chunk) for chunk in chunks]
end
