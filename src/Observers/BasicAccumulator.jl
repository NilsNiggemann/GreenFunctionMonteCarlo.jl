"""
    BasicAccumulator{T_high<:AbstractFloat} <: AbstractObserver

A basic observer struct for accumulating measurements in Monte Carlo simulations.

# Type Parameters
- `T_high<:AbstractFloat`: The floating-point type used for high-precision accumulation.

# Description
A structure that represents a basic accumulator for observing
energy and weights during simulations. Instead of storing observables at each step, the expectation values are computed on the fly, reducing the storage requirements significantly. Note that it is advisable to use a good guess of the average weight in the propagator to reduce numerical precision loss.

# Fields
- `TotalWeights`: Circular buffer of the total (normalized) walker weight at each recent step, indexed by propagation step `n`.
- `energies`: Circular buffer of the local energy estimator at each recent step.
- `Gnps`: Circular buffer holding the `Gnp[n,p]` products used to reweight energies accumulated over `p` projection steps starting from step `n`.
- `reconfigurationTable`: Circular buffer recording, for each walker and recent step, the index of the ancestor it was reconfigured from.
- `PopulationMatrix`: Circular buffer of walker population/survival information per projection order.
- `en_numerator`: Per-bin running numerator of the energy estimator, one row per projection order `p`.
- `Gnp_denominator`: Per-bin running denominator (normalization) matching `en_numerator`.
- `weight_normalization`: Scalar (0-dimensional array) normalization factor applied to `TotalWeights` to improve floating point precision.
- `bin_elements`: Number of simulation steps accumulated into a single bin before moving on to the next one.

# See Also
- [`AbstractObserver`](@ref)
"""
struct BasicAccumulator{T_high<:AbstractFloat} <: AbstractObserver
    TotalWeights::CircularArrays.CircularVector{T_high, Vector{T_high}}
    energies::CircularArrays.CircularVector{T_high, Vector{T_high}}
    Gnps::CircularArrays.CircularMatrix{T_high, Matrix{T_high}}
    reconfigurationTable::CircularArrays.CircularMatrix{Int, Matrix{Int}}
    PopulationMatrix::CircularArrays.CircularMatrix{Int, Matrix{Int}}
    en_numerator::Matrix{T_high}
    Gnp_denominator::Matrix{T_high}
    weight_normalization::Array{T_high,0}
    bin_elements::Int
end

"""
    BasicAccumulator(filename, m_proj::Integer, NWalkers::Integer; weight_normalization=1.)

Create a `BasicAccumulator` object for accumulating observables in a Green Function Monte Carlo simulation.

# Arguments
- `filename`: The name of the file where accumulated data will be stored. Providing nothing will create an in-memory accumulator.
- `m_proj::Integer`: The projection quantum number or similar parameter relevant to the simulation.
- `NWalkers::Integer`: The number of walkers used in the Monte Carlo simulation.

# Keyword Arguments
- `weight_normalization`: (default = 1.0) A normalization factor applied to the weights of the walkers. Only used to improve floating point precision.
- `num_bins`: (default = 1) The number of bins over which the accumulated energy numerator/denominator are stored separately, e.g. for error estimation via bunching.
- `bin_elements`: (default = `typemax(Int)`) The number of simulation steps accumulated into each bin before subsequent steps are assigned to the next bin.

# Returns
A `BasicAccumulator` instance configured with the specified parameters.
"""
function BasicAccumulator(filename,m_proj::Integer,NWalkers::Integer;weight_normalization=1.,num_bins=1,bin_elements=typemax(Int))
    p_proj = 2m_proj

    energies = CircularArrays.CircularArray(zeros(p_proj))
    TotalWeights = CircularArrays.CircularArray(zeros(p_proj))
    reconfigurationTable = CircularArrays.CircularArray(zeros(Int,NWalkers,p_proj))
    PopulationMatrix = CircularArrays.CircularArray(zeros(Int,NWalkers,m_proj))
    Gnps = CircularArrays.CircularArray(zeros(Float64,p_proj,p_proj))

    en_numerator = maybe_MMap_array(filename,"en_numerator",Float64,(m_proj,num_bins))
    Gnp_denominator = maybe_MMap_array(filename,"Gnp_denominator",Float64,(m_proj,num_bins))
    weight_normalization_arr = Array{Float64,0}(undef)
    weight_normalization_arr .= weight_normalization

    return BasicAccumulator(TotalWeights,energies,Gnps,reconfigurationTable,PopulationMatrix,en_numerator,Gnp_denominator,weight_normalization_arr,bin_elements)
end
BasicAccumulator(m_proj::Integer,NWalkers::Integer;kwargs...) = BasicAccumulator(nothing,m_proj,NWalkers;kwargs...)

"""
    set_zero!(A::AbstractArray)

Fill the array `A` in place with zeros of its element type.
"""
set_zero!(A::AbstractArray{T}) where T = fill!(A,zero(T))

"""
    get_num_bins(Observables::BasicAccumulator)

Return the number of bins used by `Observables` to accumulate the energy numerator and denominator.
"""
get_num_bins(Observables::BasicAccumulator) = size(Observables.en_numerator,2)

"""
    get_bin_elements(Observables::BasicAccumulator)

Return the number of simulation steps accumulated into each bin of `Observables`.
"""
get_bin_elements(Observables::BasicAccumulator) = Observables.bin_elements

"""
    projection_order(Observables::BasicAccumulator)

Return the maximal projection order `m_proj` that `Observables` was configured with.
"""
projection_order(Observables::BasicAccumulator) = size(Observables.Gnps,2) ÷2

"""
    reset_accumulator!(Observables::BasicAccumulator; hard_reset=true)

Reset the circular buffers of `Observables` (weights, energies, `Gnps`, reconfiguration and population tables) to zero.

If `hard_reset=true` (the default), also zero out the accumulated `en_numerator` and `Gnp_denominator`, discarding all previously accumulated statistics. Pass `hard_reset=false` to only reset the per-step circular buffers while keeping the running energy accumulation intact.

Returns `Observables`.
"""
function reset_accumulator!(Observables::BasicAccumulator;hard_reset=true)
    set_zero!(Observables.TotalWeights)
    set_zero!(Observables.energies)
    set_zero!(Observables.reconfigurationTable)
    set_zero!(Observables.PopulationMatrix)
    set_zero!(Observables.Gnps)
    if hard_reset 
        set_zero!(Observables.en_numerator)
        set_zero!(Observables.Gnp_denominator)
    end
    return Observables
end

"""
    get_bin_index(n, num_bins, bin_elements)
    get_bin_index(n, Observables::BasicAccumulator)

Compute the bin index that simulation step `n` (1-based) belongs to, given `bin_elements` steps per bin.

If the computed index would exceed `num_bins`, a warning is issued (once) and the last bin index is returned instead, so that later steps are folded into the final bin rather than causing an out-of-bounds access.
"""
function get_bin_index(n,num_bins,bin_elements)
    n < 1 && throw(ArgumentError("n must be greater than or equal to 1 but got n = $n"))
    n_bin_index = (n-1) ÷ bin_elements + 1
    if n_bin_index > num_bins
        @warn maxlog=1 "n_bin_index exceeds the number of bins. Using the last bin." n_bin_index num_bins
        n_bin_index = num_bins
    end
    return n_bin_index
end
get_bin_index(n,Observables::BasicAccumulator) = get_bin_index(n,get_num_bins(Observables),get_bin_elements(Observables))

"""
    saveObservables_before!(Observables::BasicAccumulator, n, Walkers, H, reconfiguration)

Update `Observables` with the measurements for propagation step `n`, to be called before the walkers are reconfigured/propagated to step `n+1`.

Updates the local energy and total weight at step `n`, propagates the `Gnp` reweighting factors, and accumulates the energy numerator and `Gnp` denominator into the appropriate bin (as determined by [`get_bin_index`](@ref)).
"""
function saveObservables_before!(Observables::BasicAccumulator,n,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Hxx = get_diagonal(H)
    energies = Observables.energies
    TotalWeights = Observables.TotalWeights

    update_energies_TotalWeights!(energies,TotalWeights,n,Walkers,Hxx)
    TotalWeights[n] /= Observables.weight_normalization[]

    Gnps = Observables.Gnps
    en_numerator = Observables.en_numerator
    Gnp_denominator = Observables.Gnp_denominator

    updateGnp!(Gnps,TotalWeights,n)
    NSites = length(getConfig(Walkers,1))
    bin_index = get_bin_index(n,Observables)
    @views getEnergy_step!(en_numerator[:,bin_index],Gnp_denominator[:,bin_index],Gnps,energies,n,NSites)
    return nothing
end

"""
    getEnergy_step!(en_numerator, Gnp_denominator, Gnp, localEnergies, n, NSites)

Accumulate the contribution of step `n` into `en_numerator` and `Gnp_denominator`, for every projection order `p < n`.

For each valid `p`, adds `Gnp[n,p] * localEnergies[n] / NSites` to `en_numerator[p]` and `Gnp[n,p]` to `Gnp_denominator[p]`, so that `en_numerator ./ Gnp_denominator` estimates the (mixed) energy per site at projection order `p`.
"""
function getEnergy_step!(en_numerator::AbstractVector,Gnp_denominator::AbstractVector,Gnp::CircularArrays.CircularMatrix,localEnergies::AbstractVector,n::Integer,NSites::Integer)
    Nsites⁻¹ = 1/NSites
    for p in eachindex(en_numerator)
        n > p || continue
        en_numerator[p] += Gnp[n,p]*localEnergies[n]*Nsites⁻¹
        Gnp_denominator[p] += Gnp[n,p]
    end
    return en_numerator
end

"""
    updateGnp!(Gnp, TotalWeights, n)

Update column `n` of the circular matrix `Gnp` with the reweighting factors accumulated over `p` projection steps ending at step `n`.

`Gnp[n,1]` records whether the weight at step `n` is positive, `Gnp[n,2]` stores the total weight at step `n`, and for `p >= 3`, `Gnp[n,p] = Gnp[n-1,p-1] * Gnp[n,2]` recursively accumulates the product of weights over the last `p-1` steps.
"""
function updateGnp!(Gnp,TotalWeights,n)
    nMax,pMax = size(Gnp)
    Gnp[n,1] = TotalWeights[n]>0 #zero projection order
    Gnp[n,2] = TotalWeights[n]

    for p in 3:pMax
        if n < p
            # Gnp[n,p] = 0
            continue
        end
        Gnp[n,p] = Gnp[n-1,p-1]*Gnp[n,2]
    end
    return
end

"""
    saveObservables_after!(Observables::BasicAccumulator, i, Walkers, H, reconfiguration)

Update `Observables` after walkers have been reconfigured at step `i`, storing the resulting reconfiguration (ancestry) list in `Observables.reconfigurationTable[:,i]`.
"""
function saveObservables_after!(Observables::BasicAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Observables.reconfigurationTable[:,i] .= get_reconfigurationList(reconfiguration)
    return nothing
end


"""
    get_energy_from_accumulator(en_numerator::AbstractMatrix, Gnp_denominator::AbstractMatrix, bin_indices::AbstractVector)

Compute the energy estimate (at each projection order) obtained from the bins in `bin_indices`, given the raw `en_numerator`/`Gnp_denominator` arrays as accumulated by a [`BasicAccumulator`](@ref) (each shaped `(m_proj, num_bins)`).

Averages the two arrays over the given bins (normalized by the mean zero-order denominator to improve numerical precision) and returns their ratio as a vector over projection orders.
"""
@views function get_energy_from_accumulator(en_numerator::AbstractMatrix,Gnp_denominator::AbstractMatrix,bin_indices::AbstractVector)
    Normalization = Statistics.mean(Gnp_denominator[1,:])

    En_num = zeros(eltype(en_numerator),size(en_numerator,1))
    Gnp_denom = zeros(eltype(Gnp_denominator),size(Gnp_denominator,1))

    for bin_idx in bin_indices
        En_num .+= en_numerator[:,bin_idx] ./ Normalization
        Gnp_denom .+= Gnp_denominator[:,bin_idx] ./ Normalization
    end
    En_num ./= Gnp_denom
    return En_num
end
"""
    get_energy_from_accumulator(en_numerator::AbstractMatrix, Gnp_denominator::AbstractMatrix)

Compute the energy estimate for each individual bin separately (i.e. without bunching), returning a vector of per-bin energy vectors.
"""
get_energy_from_accumulator(en_numerator::AbstractMatrix,Gnp_denominator::AbstractMatrix) = [get_energy_from_accumulator(en_numerator,Gnp_denominator,idx:idx) for idx in axes(Gnp_denominator,2)]

"""
    get_energy_from_accumulator_bunching(en_numerator::AbstractMatrix, Gnp_denominator::AbstractMatrix, n_bunch::Integer; kwargs...)

Compute the energy by bunching together `n_bunch` groups of bins from the raw `en_numerator`/`Gnp_denominator` arrays.

# Arguments
- `en_numerator`, `Gnp_denominator`: the accumulated numerator/denominator arrays, each shaped `(m_proj, num_bins)`, as accumulated by a [`BasicAccumulator`](@ref).
- `n_bunch::Integer`: the number of bunches to divide the bins into. For no bunching, pass `n_bunch=1`.
- `kwargs...`: additional keyword arguments forwarded to `ChunkSplitters.chunks`, such as `size` or `split`.

# Returns
A vector with one energy estimate (a vector over projection orders) per bunch.
"""
function get_energy_from_accumulator_bunching(en_numerator::AbstractMatrix,Gnp_denominator::AbstractMatrix,n_bunch::Integer;kwargs...)
    chunks = ChunkSplitters.chunks(axes(Gnp_denominator,2), size = n_bunch, split = ChunkSplitters.Consecutive();kwargs...)
    return [
        get_energy_from_accumulator(en_numerator,Gnp_denominator,chunk)
        for chunk in chunks
    ]
end

"""
    get_energy_from_accumulator(Obs::Union{BasicAccumulator,NamedTuple}, bin_indices::AbstractVector)
    get_energy_from_accumulator(Obs::Union{BasicAccumulator,NamedTuple})
    get_energy_from_accumulator_bunching(Obs::Union{BasicAccumulator,NamedTuple}, n_bunch::Integer; kwargs...)

Extract `en_numerator`/`Gnp_denominator` from `Obs` — a live [`BasicAccumulator`](@ref), or a `NamedTuple` with fields of the same name (e.g. built from an HDF5 file's `en_numerator`/`Gnp_denominator` datasets) — and forward to the array-based methods above.
"""
get_energy_from_accumulator(Obs::Union{BasicAccumulator,NamedTuple},bin_indices::AbstractVector) = get_energy_from_accumulator(Obs.en_numerator,Obs.Gnp_denominator,bin_indices)
get_energy_from_accumulator(Obs::Union{BasicAccumulator,NamedTuple}) = get_energy_from_accumulator(Obs.en_numerator,Obs.Gnp_denominator)
get_energy_from_accumulator_bunching(Obs::Union{BasicAccumulator,NamedTuple},n_bunch::Integer;kwargs...) = get_energy_from_accumulator_bunching(Obs.en_numerator,Obs.Gnp_denominator,n_bunch;kwargs...)

"""
    get_energy_from_accumulator_bunching(filename::AbstractString, n_bunch::Integer; kwargs...)

Read the `en_numerator`/`Gnp_denominator` datasets directly from the HDF5 file written by a [`BasicAccumulator`](@ref) (i.e. constructed as `BasicAccumulator(filename, ...)`), and bunch them into `n_bunch` groups.
"""
function get_energy_from_accumulator_bunching(filename::AbstractString,n_bunch::Integer;kwargs...)
    en_numerator,Gnp_denominator = HDF5.h5open(filename,"r") do file
        read(file["en_numerator"]),read(file["Gnp_denominator"])
    end
    return get_energy_from_accumulator_bunching(en_numerator,Gnp_denominator,n_bunch;kwargs...)
end

"""
    log_observable(O::BasicAccumulator, i)

Return the tuple of log entries reported for step `i`: the walker survival ratio, the local energy, and the average weight.

Used to populate progress/log output during a simulation run.
"""
log_observable(O::BasicAccumulator,i) = (log_walker_survival_ratio(O.reconfigurationTable,i),log_Obs_energy(O,i),log_Obs_weights(O,i))

"""
    log_Obs_energy(O::BasicAccumulator, i)

Return a `("eloc", value)` pair with the local energy at step `i`, formatted for logging.
"""
log_Obs_energy(O::BasicAccumulator,i) = ("eloc",strd(O.energies[i]))

"""
    log_Obs_weights(O::BasicAccumulator, i)

Return a `("w_avg", value)` pair with the total walker weight at step `i`, formatted for logging.
"""
log_Obs_weights(O::BasicAccumulator,i) = ("w_avg",strd(O.TotalWeights[i]))
