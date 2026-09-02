"""
    SymmetrizedSpinCorrelations(inds_Rij, t = Float32; normalize = true) -> SymmetrizedSpinCorrelations{t}

Constructs a symmetry-reduced version of [`SpinCorrelations`](@ref): instead of measuring every
pair ``\\langle S^z_i S^z_j\\rangle`` separately, the correlations are averaged over groups of
symmetry-equivalent site pairs,
```math
C(R) = \\frac{1}{N_R} \\sum_{(i,j) \\in R} S^z_i S^z_j ,
```
where each group `R` (e.g. all pairs with the same distance vector ``r_i - r_j = R``) is supplied
by the user. The symmetry reduction itself is deliberately left outside of this package, so that
partially broken symmetries, sublattice resolution, etc. can be handled by the caller.

Since the buffer holds one entry per group instead of ``N_\\mathrm{sites}(N_\\mathrm{sites}+1)/2``
entries, this is considerably cheaper to accumulate than [`SpinCorrelations`](@ref) — the numerator
arrays of [`ObservableAccumulator`](@ref)/[`LazyObservableAccumulator`](@ref) scale linearly in the
number of measured values.

# Arguments
- `inds_Rij`: An iterable of groups, where each group is an iterable of site pairs `(i,j)` that are
  summed and averaged into a single output value. Pairs may be given as `Tuple`s, `Pair`s,
  `CartesianIndex{2}` or two-element vectors, e.g.
  `[[(1,1),(2,2),(3,3)], [(1,2),(2,3),(3,1)]]`.
- `t`: (optional) The numeric type of the internal buffer, defaulting to `Float32`. Note that,
  unlike for [`SpinCorrelations`](@ref), an integer type is only meaningful together with
  `normalize = false`, since the averaged correlations are not integers.
- `normalize`: (optional, keyword) If `true` (the default) each group is divided by its number of
  pairs ``N_R``; if `false`, the plain sum over the group is measured.

# Returns
A `SymmetrizedSpinCorrelations` object whose buffer has one entry per group of `inds_Rij`, in the
order in which the groups were given.

# Example
```julia
L = 8 # periodic chain: group the pairs by their distance R along the chain
inds_Rij = [[(i, mod1(i + R, L)) for i in 1:L] for R in 0:L-1]
Observable = SymmetrizedSpinCorrelations(inds_Rij)
```

# See also
- [`SpinCorrelations`](@ref)
- [`LazyObservableAccumulator`](@ref)
- [`ObservableAccumulator`](@ref)
"""
struct SymmetrizedSpinCorrelations{T<:Real} <: AbstractObservable
    obsBuffer::Vector{T}
    i_inds::Vector{Int}
    j_inds::Vector{Int}
    offsets::Vector{Int}
    norms::Vector{T}
end

_pair_ij(ij::Tuple{Integer,Integer}) = (Int(ij[1]), Int(ij[2]))
_pair_ij(ij::Pair{<:Integer,<:Integer}) = (Int(first(ij)), Int(last(ij)))
_pair_ij(ij::CartesianIndex{2}) = (Int(ij[1]), Int(ij[2]))
function _pair_ij(ij::AbstractVector{<:Integer})
    length(ij) == 2 || throw(ArgumentError("each site pair must consist of exactly two indices, got $(length(ij))"))
    return (Int(ij[begin]), Int(ij[begin+1]))
end
_pair_ij(ij) = throw(ArgumentError("cannot interpret $(ij)::$(typeof(ij)) as a site pair (i,j). Use a Tuple, Pair, CartesianIndex{2} or a two-element vector."))

function SymmetrizedSpinCorrelations(inds_Rij, t = Float32; normalize = true)
    if normalize && !(t <: AbstractFloat)
        throw(ArgumentError("normalize = true requires a floating point buffer type, got $t. Pass normalize = false to measure the unnormalized sum over each group instead."))
    end
    groups = collect(inds_Rij)
    isempty(groups) && throw(ArgumentError("inds_Rij must contain at least one group of site pairs"))

    NGroups = length(groups)
    offsets = Vector{Int}(undef, NGroups + 1)
    offsets[begin] = 1
    norms = Vector{t}(undef, NGroups)
    i_inds = Int[]
    j_inds = Int[]

    for (k, group) in enumerate(groups)
        ij_pairs = collect(group)
        isempty(ij_pairs) && throw(ArgumentError("group $k of inds_Rij is empty: every group must contain at least one site pair (i,j)"))
        for ij in ij_pairs
            (i, j) = _pair_ij(ij)
            (i > 0 && j > 0) || throw(ArgumentError("site indices must be positive, got ($i,$j) in group $k of inds_Rij"))
            push!(i_inds, i)
            push!(j_inds, j)
        end
        offsets[k+1] = offsets[k] + length(ij_pairs)
        norms[k] = normalize ? one(t) / length(ij_pairs) : one(t)
    end

    return SymmetrizedSpinCorrelations(zeros(t, NGroups), i_inds, j_inds, offsets, norms)
end

Base.copy(O::SymmetrizedSpinCorrelations) = SymmetrizedSpinCorrelations(copy(O.obsBuffer), copy(O.i_inds), copy(O.j_inds), copy(O.offsets), copy(O.norms))
obs(O::SymmetrizedSpinCorrelations) = O.obsBuffer

function (O::SymmetrizedSpinCorrelations)(out, config::AbstractVector)
    for k in eachindex(O.norms)
        # the pairs of group k occupy the contiguous range p_start:p_stop of i_inds/j_inds
        p_start = O.offsets[k]
        p_stop = O.offsets[k+1] - 1
        # accumulate in Float64 independently of the buffer type: a group may contain many pairs,
        # and the sum is divided by their number before it is written to the buffer.
        Corr = 0.

        LoopVectorization.@turbo for p in p_start:p_stop
            i = O.i_inds[p]
            j = O.j_inds[p]
            xi = config[i]
            xj = config[j]
            Corr += (1. - 2. * xi) * (1. - 2. * xj)
        end

        out[k] = Corr * O.norms[k]
    end
    return out
end

(O::SymmetrizedSpinCorrelations)(out, config::BosonConfig) = O(out, parent(config))

function average_obs_walkers(O::SymmetrizedSpinCorrelations, obs_idx::Integer, walker_confs::AbstractMatrix, WalkerPopulations::AbstractVector{<:Integer})
    # accumulate in Float64: both the number of walkers and the number of pairs per group can be
    # large, and the result is added to the Float64 numerators of the accumulator anyway.
    Obs_num_i_m_b = 0.

    p_start = O.offsets[obs_idx]
    p_stop = O.offsets[obs_idx+1] - 1

    for p in p_start:p_stop
        i = O.i_inds[p]
        j = O.j_inds[p]
        Corr_ij = 0.

        LoopVectorization.@turbo for α in axes(walker_confs, 1)
            mult = WalkerPopulations[α]
            x_i = walker_confs[α, i]
            x_j = walker_confs[α, j]

            Corr_ij += (1. - 2. * x_i) * (1. - 2. * x_j) * mult
        end

        Obs_num_i_m_b += Corr_ij
    end
    return Obs_num_i_m_b * O.norms[obs_idx]
end
