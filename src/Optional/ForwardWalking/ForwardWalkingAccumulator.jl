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
existing diagonal `Gnp[n,p]` energy estimator, realized here by literal forward continuation rather
than backward reweighting. Check convergence of [`get_observable_from_accumulator`](@ref)'s output
vs. projection order `p` to validate that `mProj` is large enough.

# Requirement
`FW.BasicAcc` must be the *same* [`BasicAccumulator`](@ref) instance also included, separately, as
its own observer in the enclosing [`CombinedObserver`](@ref) — exactly as already required for
[`ObservableAccumulator`](@ref) — so that its `Gnps` table is kept up to date; `ForwardWalkingAccumulator`
only reads it, it never updates it itself.
"""
struct ForwardWalkingAccumulator{ObsType<:AbstractOffdiagonalObservable,SlotType,Prop<:AbstractPropagator,GF<:AbstractGuidingFunction,HS<:AbstractHilbertSpace,PS<:AbstractParallelizationScheme,T_high<:AbstractFloat} <: AbstractObserver
    BasicAcc::BasicAccumulator{T_high}
    Observable::ObsType
    slots::Vector{SlotType}
    cadence::Int
    mProj::Int
    propagator::Prop
    logψ::GF
    Hilbert::HS
    parallelization::PS
    w_avg_estimate::Float64
    weights_buf::Vector{Float64}
    rng::Random.AbstractRNG
    Num::Matrix{Float64}
    Denom::Matrix{Float64}
end

"""
    ForwardWalkingAccumulator(filename, Observable::AbstractOffdiagonalObservable, BasicAcc::BasicAccumulator,
                               conf, logψ::AbstractGuidingFunction, H::AbstractSignFreeOperator, Hilbert::AbstractHilbertSpace,
                               propagator::AbstractPropagator, mProj::Integer, cadence::Integer, Nsub::Integer;
                               w_avg_estimate=1., parallelization=SingleThreaded(), rng=Random.default_rng(), num_bins=get_num_bins(BasicAcc))

Construct a [`ForwardWalkingAccumulator`](@ref) for measuring `⟨Observable⟩` via streaming forward walking.

# Arguments
- `filename`: HDF5 file to persist `Num`/`Denom` to (mirroring `BasicAccumulator`/`ObservableAccumulator`), or `nothing` for an in-memory-only accumulator.
- `Observable::AbstractOffdiagonalObservable`: the off-diagonal operator to measure.
- `BasicAcc::BasicAccumulator`: the accumulator supplying the backward `Gnp` denominator. Must also be included, separately, in the enclosing `CombinedObserver` (see the type docstring).
- `conf`: an exemplary configuration, used to allocate each forward-walking sub-population.
- `logψ`, `H`, `Hilbert`, `propagator`: the same guiding function, Hamiltonian, Hilbert space, and propagator used by the main simulation — sub-populations are continued with *ordinary* dynamics after being seeded.
- `mProj::Integer`: number of forward-walking continuation steps (the projection-order sweep goes from `p=0` to `p=mProj`). Must satisfy `mProj < 2*m_proj_Basic`, where `m_proj_Basic` is `BasicAcc`'s own `m_proj` (asserted at construction).
- `cadence::Integer`: seed a new forward-walking sub-population every `cadence` main-loop steps.
- `Nsub::Integer`: size of each forward-walking sub-population. Must equal the main population's `NWalkers` (v1: 1:1 seeding, no resampling).

# Keyword Arguments
- `w_avg_estimate`: average branching-weight estimate used to rescale sub-population weights (same role as for `ContinuousTimePropagator`/`DiscretePropagator`). Default `1.0`.
- `parallelization`: parallelization scheme for advancing sub-populations. Default `SingleThreaded()`.
- `rng`: random number generator used for seeding/advancing sub-populations. Note that `saveObservables_before!`/`_after!` do not receive the main loop's RNG, so this accumulator owns its own stream.
- `num_bins`: number of bins for error-bar bunching (default: `BasicAcc`'s own).
"""
function ForwardWalkingAccumulator(filename, Observable::AbstractOffdiagonalObservable, BasicAcc::BasicAccumulator,
                                    conf, logψ::AbstractGuidingFunction, H::AbstractSignFreeOperator, Hilbert::AbstractHilbertSpace,
                                    propagator::AbstractPropagator, mProj::Integer, cadence::Integer, Nsub::Integer;
                                    w_avg_estimate::Real = 1.,
                                    parallelization::AbstractParallelizationScheme = SingleThreaded(),
                                    rng::Random.AbstractRNG = Random.default_rng(),
                                    num_bins::Integer = get_num_bins(BasicAcc))

    m_proj_Basic = size(BasicAcc.Gnps, 1) ÷ 2
    @assert mProj < 2 * m_proj_Basic "ForwardWalkingAccumulator requires mProj < 2*m_proj_Basic (the wrapped BasicAccumulator's own m_proj), got mProj=$mProj, 2*m_proj_Basic=$(2*m_proj_Basic); otherwise Gnps's circular buffer wraps before a forward-walking slot finishes reading it, silently corrupting results."

    n_slots = fld(mProj, cadence) + 1
    slots = [allocate_forward_walking_slot(conf, logψ, Nsub, H) for _ in 1:n_slots]

    weights_buf = zeros(n_moves(Observable))

    Num = maybe_MMap_array(filename, "$(_type_stripped(Observable))_ForwardWalking_numerator", Float64, (mProj + 1, num_bins))
    Denom = maybe_MMap_array(filename, "$(_type_stripped(Observable))_ForwardWalking_denominator", Float64, (mProj + 1, num_bins))

    return ForwardWalkingAccumulator(BasicAcc, Observable, slots, Int(cadence), Int(mProj), propagator, logψ, Hilbert,
                                      parallelization, Float64(w_avg_estimate), weights_buf, rng, Num, Denom)
end

function ForwardWalkingAccumulator(Observable::AbstractOffdiagonalObservable, BasicAcc::BasicAccumulator, args...; kwargs...)
    return ForwardWalkingAccumulator(nothing, Observable, BasicAcc, args...; kwargs...)
end

function saveObservables_before!(FW::ForwardWalkingAccumulator, n, Walkers::AbstractWalkerEnsemble, H::AbstractSignFreeOperator, reconfiguration::AbstractReconfigurationScheme)
    bin_index = get_bin_index(n, FW.BasicAcc)

    for slot in FW.slots
        slot.active || continue
        advance_slot!(slot, FW.propagator, FW.logψ, H, FW.Hilbert, FW.w_avg_estimate, FW.parallelization, FW.rng)
        c = slot.step + 1
        # Mirrors getEnergy_step!'s `n > p` guard (BasicAccumulator.jl): Gnps[n_seed,c] is only a
        # meaningful backward-projected weight once n_seed has at least c steps of history behind it;
        # for seeds occurring in the first ~mProj main-loop steps this skips the not-yet-valid columns
        # (negligible for NSteps >> mProj, exactly as for the diagonal energy estimator's own warm-up).
        if slot.n_seed > c
            FW.Num[c, bin_index] += slot.seed_weight * slot.cum_weight
            FW.Denom[c, bin_index] += FW.BasicAcc.Gnps[slot.n_seed, c]
        end
        slot.step == FW.mProj && (slot.active = false)
    end

    if mod(n, FW.cadence) == 0
        slot = _find_free_slot(FW.slots)
        seed_from_main_population!(slot, Walkers, FW.Observable, FW.logψ, FW.Hilbert, n, FW.weights_buf, FW.rng)
        if n > 1
            FW.Num[1, bin_index] += slot.seed_weight
            FW.Denom[1, bin_index] += FW.BasicAcc.Gnps[n, 1]
        end
    end

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

"""
    get_observable_from_accumulator_bunching(FW::ForwardWalkingAccumulator, n_bunch::Integer; kwargs...)

Like [`get_observable_from_accumulator`](@ref), but groups the accumulated bins into `n_bunch` groups
first (for statistical error analysis), mirroring [`get_energy_from_accumulator_bunching`](@ref).
"""
function get_observable_from_accumulator_bunching(FW::ForwardWalkingAccumulator, n_bunch::Integer; kwargs...)
    chunks = ChunkSplitters.chunks(axes(FW.Denom, 2), size = n_bunch, split = ChunkSplitters.Consecutive(); kwargs...)
    return [get_observable_from_accumulator(FW, chunk) for chunk in chunks]
end
