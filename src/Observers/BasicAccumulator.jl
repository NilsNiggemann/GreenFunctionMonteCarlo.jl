"""
    BasicAccumulator{T_high<:AbstractFloat} <: AbstractObserver

A basic observer struct for accumulating measurements in Monte Carlo simulations.

# Type Parameters
- `T_high<:AbstractFloat`: The floating-point type used for high-precision accumulation.

# Description
A structure that represents a basic accumulator for observing
energy and weights during simulations. Instead of storing observables at each step, the expectation values are computed on the fly, reducing the storage requirements significantly. Note that it is advisable to use a good guess of the average weight in the propagator to reduce numerical precision loss.

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
set_zero!(A::AbstractArray{T}) where T = fill!(A,zero(T))
get_num_bins(Observables::BasicAccumulator) = size(Observables.en_numerator,2)
get_bin_elements(Observables::BasicAccumulator) = Observables.bin_elements

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

BasicAccumulator(m_proj::Integer,NWalkers::Integer;kwargs...) = BasicAccumulator(nothing,m_proj,NWalkers;kwargs...)

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

function getEnergy_step!(en_numerator::AbstractVector,Gnp_denominator::AbstractVector,Gnp::CircularArrays.CircularMatrix,localEnergies::AbstractVector,n::Integer,NSites::Integer)
    Nsites⁻¹ = 1/NSites
    for p in eachindex(en_numerator)
        n > p || continue
        en_numerator[p] += Gnp[n,p]*localEnergies[n]*Nsites⁻¹
        Gnp_denominator[p] += Gnp[n,p]
    end
    return en_numerator
end

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

function saveObservables_after!(Observables::BasicAccumulator,i,Walkers::AbstractWalkerEnsemble,H::AbstractSignFreeOperator,reconfiguration::AbstractReconfigurationScheme)
    Observables.reconfigurationTable[:,i] .= get_reconfigurationList(reconfiguration)
    return nothing
end


"""
    get_energy_from_accumulator_bunching(Observables::BasicAccumulator, n_bunch::Integer; kwargs...)

Compute the energy by bunching together observable accumulators.

# Arguments
- `Observables::BasicAccumulator`: The accumulator containing observable measurements.
- `n_bunch::Integer`: The number of bunches to divide the data into for statistical analysis. For no bunching, pass `n_bunch=1`.
- `kwargs...`: Additional keyword arguments for customization of the chunking process, such as `size` or `split`.

# Returns
- The computed energy, possibly with statistical error estimates depending on implementation.

# Description
This function processes the accumulated observables by grouping (bunching) the data into `n_bunch` groups. It then calculates the energy, which can be used for error analysis or to reduce autocorrelation effects in Monte Carlo simulations.
"""
function get_energy_from_accumulator_bunching(Observables::Union{BasicAccumulator,<:NamedTuple},n_bunch::Integer;kwargs...)
    chunks = ChunkSplitters.chunks(axes(Observables.Gnp_denominator,2), size = n_bunch, split = ChunkSplitters.Consecutive();kwargs...)
    return [
        get_energy_from_accumulator(Observables,chunk)
        for chunk in chunks
    ]
end
    
@views function get_energy_from_accumulator(Obs::Union{BasicAccumulator,<:NamedTuple},bin_indices::AbstractVector)
    Normalization = Statistics.mean(Obs.Gnp_denominator[1,:])

    En_num = zeros(eltype(Obs.en_numerator),size(Obs.en_numerator,1))
    Gnp_denom = zeros(eltype(Obs.Gnp_denominator),size(Obs.Gnp_denominator,1))

    for bin_idx in bin_indices
        En_num .+= Obs.en_numerator[:,bin_idx] ./ Normalization
        Gnp_denom .+= Obs.Gnp_denominator[:,bin_idx] ./ Normalization 
    end
    En_num ./= Gnp_denom
    return En_num
end
get_energy_from_accumulator(Observables) = [get_energy_from_accumulator(Observables,idx:idx) for idx in axes(Observables.Obs_denominators,2)]

log_observable(O::BasicAccumulator,i) = (log_walker_survival_ratio(O.reconfigurationTable,i),log_Obs_energy(O,i),log_Obs_weights(O,i))
log_Obs_energy(O::BasicAccumulator,i) = ("eloc",strd(O.energies[i]))
log_Obs_weights(O::BasicAccumulator,i) = ("w_avg",strd(O.TotalWeights[i]))
