@doc raw"""
    ClusterCPS{T<:AbstractFloat,M} <: AbstractGuidingFunction

Cluster correlator product state (CPS), also known as entangled-plaquette state: a product of
correlators, one per cluster of sites, each tabulated over all local states of its cluster,
```math
\log\psi(x) = \sum_c T_{t(c)}\big[I_c(x)\big], \qquad
I_c(x) = \sum_{k=1}^{L_c} \big(x_{s_{c,k}} - n_\mathrm{min}\big)\, d^{\,k-1},
```
where cluster ``c = (s_{c,1}, \dots, s_{c,L_c})`` is an ordered list of ``L_c`` distinct sites,
``d`` is the local dimension (`local_dim`, 2 for hard-core bosons / spin-1/2), ``n_\mathrm{min}``
the smallest occupation (`min_occupation`), and ``T_t`` a table with ``d^{L_c}`` entries
(1-based Julia index ``I_c + 1``). The first site of a cluster is the least significant digit.

Clusters may have different sizes, may overlap, and need not cover all sites. Several clusters
may share one table (e.g. all translates of a plaquette), which ties their parameters together.
Table entries may be `-Inf` to forbid a cluster pattern (``\psi = 0``); walkers must start
in configurations with finite ``\log\psi``.

# Type Parameters
- `T<:AbstractFloat`: element type of the tables (integer tables are converted with `float`).
- `M`: element type of the precompiled `moves` (`Nothing` if none were given).

# Constructors
- `ClusterCPS(clusters, tables, cluster_table; local_dim = 2, min_occupation = 0, moves = nothing)`:
  cluster `c` uses `tables[cluster_table[c]]`.
- `ClusterCPS(clusters, table; kwargs...)`: all clusters share `table` (all must have equal size).

`clusters` is a vector of site-index vectors, or a matrix whose columns are the clusters.

# Keyword Arguments
- `local_dim::Integer = 2`: number of local states ``d``; every site of a cluster must have an
  occupation in `min_occupation:min_occupation+local_dim-1` (checked, also for the target of a move).
- `min_occupation::Integer = 0`: occupation that maps to digit 0.
- `moves = nothing`: optionally the moves of the Hamiltonian (e.g. `H.moves`). Each of them is then
  precompiled into the list of clusters it touches and their pattern shifts, which makes
  `log_psi_diff` about 1.5-2× faster. Supported are moves with a state-independent `move_dx`
  (e.g. `SparseMove`) and `FlipMove`s (requires `local_dim = 2`). Moves not in this list, or not
  identical (`===`) to one of them, are still handled exactly by the generic code path.

# Cost
`log_psi_diff` touches only the clusters containing an affected site of the move, i.e.
``\mathcal O(\sum_{i \in \text{move}} \#\{c \ni i\})`` operations, independent of the system size.
The buffer stores the integer patterns ``I_c(x)``, which are updated exactly (no round-off drift).
The buffer also holds scratch space used by the generic path of `log_psi_diff`, so one buffer must
not be used by several threads at once (as for every walker buffer in this package).

# Example
```julia
N = 8
Hilbert = BosonHilbertSpace(N, HardCoreConstraint())
bonds = [[i, mod1(i + 1, N)] for i in 1:N]          # nearest-neighbour pairs, shared table
logψ = ClusterCPS(bonds, [0.0, 0.1, 0.1, -0.3])    # patterns (0,0), (1,0), (0,1), (1,1)
```

See also [`Jastrow`](@ref), [`get_params`](@ref), [`GreenFunctionMonteCarlo.get_tables`](@ref).
"""
struct ClusterCPS{T<:AbstractFloat,M} <: AbstractGuidingFunction
    params::Vector{T}                       # all tables, concatenated
    table_ranges::Vector{UnitRange{Int}}    # table t -> its entries in params
    cluster_table::Vector{Int}              # cluster c -> table t
    cluster_offset::Vector{Int}             # params index of pattern 0 of cluster c
    cluster_ptr::Vector{Int}                # sites of cluster c: cluster_sites[cluster_ptr[c]:cluster_ptr[c+1]-1]
    cluster_sites::Vector{Int}
    site_ptr::Vector{Int}                   # clusters containing site i: site_cluster[site_ptr[i]:site_ptr[i+1]-1]
    site_cluster::Vector{Int}
    site_stride::Vector{Int}                # d^(k-1) for the slot k that site i occupies in that cluster
    local_dim::Int
    min_occupation::Int
    moves::Vector{M}                        # precompiled moves (empty if none were given)
    move_ptr::Vector{Int}                   # moves whose first affected site is i: move_ids[move_ptr[i]:move_ptr[i+1]-1]
    move_ids::Vector{Int}
    row_ptr::Vector{Int}                    # clusters touched by precompiled move j: rows row_ptr[j]:row_ptr[j+1]-1
    row_cluster::Vector{Int}
    row_shift::Vector{Int}                  # pattern shift (added, or XOR mask for a FlipMove) of that cluster
end

function ClusterCPS(clusters::AbstractVector{<:AbstractVector{<:Integer}}, tables::AbstractVector{<:AbstractVector{T}}, cluster_table::AbstractVector{<:Integer}; local_dim::Integer = 2, min_occupation::Integer = 0, moves = nothing) where {T<:Real}
    local_dim >= 2 || throw(ArgumentError("local_dim must be at least 2, got $local_dim"))
    length(cluster_table) == length(clusters) || throw(ArgumentError("cluster_table has $(length(cluster_table)) entries for $(length(clusters)) clusters"))
    for (c, cluster) in enumerate(clusters)
        isempty(cluster) && throw(ArgumentError("cluster $c is empty"))
        all(>=(1), cluster) || throw(ArgumentError("cluster $c contains a site index < 1: $cluster"))
        allunique(cluster) || throw(ArgumentError("cluster $c contains a site more than once: $cluster"))
        t = cluster_table[c]
        1 <= t <= length(tables) || throw(ArgumentError("cluster $c refers to table $t, but there are $(length(tables)) tables"))
        length(tables[t]) == big(local_dim)^length(cluster) || throw(ArgumentError("table $t has $(length(tables[t])) entries, but cluster $c of size $(length(cluster)) needs local_dim^$(length(cluster)) = $(big(local_dim)^length(cluster))"))
    end
    for (t, table) in enumerate(tables)
        all(<(Inf), table) || throw(ArgumentError("table $t contains NaN or +Inf (only -Inf, i.e. ψ = 0, is allowed)"))
    end

    params = Vector{float(T)}(reduce(vcat, tables; init = T[]))
    table_ranges = Vector{UnitRange{Int}}(undef, length(tables))
    first_entry = 1
    for (t, table) in enumerate(tables)
        table_ranges[t] = first_entry:first_entry+length(table)-1
        first_entry += length(table)
    end
    cluster_offset = [first(table_ranges[t]) for t in cluster_table]
    cluster_ptr = cumsum([1; length.(clusters)])
    cluster_sites = reduce(vcat, clusters; init = Int[])

    nsites = maximum(cluster_sites; init = 0)
    counts = zeros(Int, nsites)
    foreach(i -> counts[i] += 1, cluster_sites)
    site_ptr = cumsum([1; counts])
    site_cluster, site_stride = similar(cluster_sites), similar(cluster_sites)
    next = site_ptr[1:end-1]
    for (c, cluster) in enumerate(clusters), (k, i) in enumerate(cluster)
        site_cluster[next[i]], site_stride[next[i]] = c, local_dim^(k - 1)
        next[i] += 1
    end

    compiled = _compile_moves(isnothing(moves) ? Nothing[] : collect(moves), site_ptr, site_cluster, site_stride, local_dim, min_occupation)
    return ClusterCPS(params, table_ranges, Vector{Int}(cluster_table), cluster_offset, cluster_ptr, cluster_sites, site_ptr, site_cluster, site_stride, Int(local_dim), Int(min_occupation), compiled...)
end

ClusterCPS(clusters::AbstractVector{<:AbstractVector{<:Integer}}, table::AbstractVector{<:Real}; kwargs...) = ClusterCPS(clusters, [table], ones(Int, length(clusters)); kwargs...)
ClusterCPS(clusters::AbstractMatrix{<:Integer}, args...; kwargs...) = ClusterCPS(collect(eachcol(clusters)), args...; kwargs...)

# Coefficient of slot stride for site `a` of a precompiled move, and how the shift acts on a pattern.
_shift_coefficient(move::FlipMove, a) = 1
_shift_coefficient(move::AbstractMove, a) = Int(move_dx(move)[a])
@inline _shifted(::FlipMove, pattern, shift) = pattern ⊻ shift
@inline _shifted(::AbstractMove, pattern, shift) = pattern + shift

function _check_compilable(move::FlipMove, local_dim, min_occupation)
    (local_dim, min_occupation) == (2, 0) || throw(ArgumentError("FlipMoves can only be precompiled for local_dim = 2 and min_occupation = 0"))
end
function _check_compilable(move::AbstractMove, local_dim, min_occupation)
    hasmethod(move_dx, Tuple{typeof(move)}) || throw(ArgumentError("moves of type $(typeof(move)) have no state-independent move_dx and cannot be precompiled; pass moves = nothing"))
end

# Touched clusters and pattern shifts of every move, plus a first-site index for looking moves up.
function _compile_moves(moves::Vector{M}, site_ptr, site_cluster, site_stride, local_dim, min_occupation) where {M}
    row_ptr, row_cluster, row_shift = [1], Int[], Int[]
    first_sites = zeros(Int, length(moves))
    for (j, move) in enumerate(moves)
        _check_compilable(move, local_dim, min_occupation)
        sites = affected_sites(move)
        shifts = Dict{Int,Int}()
        for (a, i) in enumerate(sites), q in (i < length(site_ptr) ? (site_ptr[i]:site_ptr[i+1]-1) : (1:0))
            c = site_cluster[q]
            shifts[c] = get(shifts, c, 0) + _shift_coefficient(move, a) * site_stride[q]
        end
        for (c, shift) in sort!(collect(shifts))
            iszero(shift) && continue
            push!(row_cluster, c)
            push!(row_shift, shift)
        end
        push!(row_ptr, length(row_cluster) + 1)
        isempty(sites) || (first_sites[j] = first(sites))
    end
    counts = zeros(Int, maximum(first_sites; init = 0))
    foreach(i -> i > 0 && (counts[i] += 1), first_sites)
    move_ptr = cumsum([1; counts])
    move_ids = zeros(Int, length(moves) - count(iszero, first_sites))
    next = move_ptr[1:end-1]
    for (j, i) in enumerate(first_sites)
        i > 0 || continue
        move_ids[next[i]] = j
        next[i] += 1
    end
    return moves, move_ptr, move_ids, row_ptr, row_cluster, row_shift
end
_compile_moves(moves::Vector{Nothing}, args...) = (Nothing[], [1], Int[], [1], Int[], Int[])

guidingfunc_name(::ClusterCPS) = "ClusterCPS"

"""
    get_params(logψ::ClusterCPS)

All tables concatenated into one vector (a shared table appears once). Mutating it changes `logψ`
immediately; walker buffers stay valid, since they only store the cluster patterns.
"""
get_params(logψ::ClusterCPS) = logψ.params

"""
    get_tables(logψ::ClusterCPS)

Views of the individual tables into [`get_params`](@ref)`(logψ)`, in the order they were passed.
"""
get_tables(logψ::ClusterCPS) = [view(logψ.params, r) for r in logψ.table_ranges]

num_clusters(logψ::ClusterCPS) = length(logψ.cluster_offset)

@noinline _throw_local_state(i, v, logψ) = throw(ArgumentError("site $i has occupation $v outside the range $(logψ.min_occupation):$(logψ.min_occupation + logψ.local_dim - 1) covered by the ClusterCPS (local_dim = $(logψ.local_dim), min_occupation = $(logψ.min_occupation))"))

@inline function _digit(logψ::ClusterCPS, i, v)
    digit = Int(v) - logψ.min_occupation
    0 <= digit < logψ.local_dim || _throw_local_state(i, v, logψ)
    return digit
end

# Pattern I_c(x) of cluster c (0-based).
function cluster_pattern(x::AbstractArray, logψ::ClusterCPS, c::Integer)
    I, stride = 0, 1
    for p in logψ.cluster_ptr[c]:logψ.cluster_ptr[c+1]-1
        i = logψ.cluster_sites[p]
        I += _digit(logψ, i, x[i]) * stride
        stride *= logψ.local_dim
    end
    return I
end

function (logψ::ClusterCPS{T})(x::AbstractArray) where {T}
    s = zero(T)
    for c in 1:num_clusters(logψ)
        s += logψ.params[logψ.cluster_offset[c]+cluster_pattern(x, logψ, c)]
    end
    return s
end

"""
    ClusterCPSBuffer <: AbstractGuidingFunctionBuffer

Walker buffer of a [`ClusterCPS`](@ref): the current pattern of every cluster, plus scratch space
for `log_psi_diff` (`delta`, `stamp`, `touched`, `generation`), which is never copied between buffers.
"""
struct ClusterCPSBuffer <: AbstractGuidingFunctionBuffer
    pattern::Vector{Int}                  # I_c(x)
    delta::Vector{Int}                    # scratch: pattern change of cluster c under the move being evaluated
    stamp::Vector{Int}                    # scratch: delta[c] is valid iff stamp[c] == generation[]
    touched::Vector{Int}                  # scratch: clusters touched by the move being evaluated
    generation::Base.RefValue{Int}
end

function allocate_GWF_buffer(logψ::ClusterCPS, x)
    nc = num_clusters(logψ)
    Buffer = ClusterCPSBuffer(zeros(Int, nc), zeros(Int, nc), zeros(Int, nc), zeros(Int, nc), Ref(0))
    return compute_GWF_buffer!(Buffer, logψ, x)
end

function compute_GWF_buffer!(Buffer::ClusterCPSBuffer, logψ::ClusterCPS, x)
    for c in 1:num_clusters(logψ)
        Buffer.pattern[c] = cluster_pattern(x, logψ, c)
    end
    return Buffer
end

function setBuffer!(A::ClusterCPSBuffer, B::ClusterCPSBuffer)
    A.pattern .= B.pattern
    return A
end

@inline _site_range(logψ::ClusterCPS, i) = i < length(logψ.site_ptr) ? (logψ.site_ptr[i]:logψ.site_ptr[i+1]-1) : (1:0)

@inline function log_psi_diff(x::AbstractConfig, move::AbstractMove, logψ::ClusterCPS{T}, Buffer::ClusterCPSBuffer, Hilbert::AbstractHilbertSpace) where {T}
    isapplicable(x, move, Hilbert) || return T(-Inf)
    j = _compiled_move_id(logψ, move)
    iszero(j) && return _log_psi_diff_generic(x, move, logψ, Buffer)
    _covers_hilbert_space(logψ, constraint(Hilbert)) || _check_target_states(x, move, logψ)
    Δ = zero(T)
    @inbounds for r in logψ.row_ptr[j]:logψ.row_ptr[j+1]-1
        c = logψ.row_cluster[r]
        offset, old = logψ.cluster_offset[c], Buffer.pattern[c]
        Δ += logψ.params[offset+_shifted(move, old, logψ.row_shift[r])] - logψ.params[offset+old]
    end
    return Δ
end

# Index of `move` among the precompiled moves (compared with `===`), or 0.
@inline function _compiled_move_id(logψ::ClusterCPS, move::AbstractMove)
    sites = affected_sites(move)
    isempty(sites) && return 0
    i = first(sites)
    i < length(logψ.move_ptr) || return 0
    @inbounds for p in logψ.move_ptr[i]:logψ.move_ptr[i+1]-1
        j = logψ.move_ids[p]
        logψ.moves[j] === move && return j
    end
    return 0
end

# true if every occupation allowed by the constraint is a valid local state of the CPS, so targets of applicable moves need no check
_covers_hilbert_space(logψ::ClusterCPS, ::AbstractConstraint) = false
_covers_hilbert_space(logψ::ClusterCPS, ::HardCoreConstraint) = logψ.min_occupation == 0
_covers_hilbert_space(logψ::ClusterCPS, c::OccupationNumberConstraint) = logψ.min_occupation <= c.min_occupation && c.max_occupation < logψ.min_occupation + logψ.local_dim

_check_target_states(x, move::FlipMove, logψ::ClusterCPS) = nothing   # Bool configs, local_dim = 2
@inline function _check_target_states(x, move::AbstractMove, logψ::ClusterCPS)
    dx = move_dx_before(move, x)
    for (a, i) in enumerate(affected_sites(move))
        isempty(_site_range(logψ, i)) || _digit(logψ, i, x[i] + dx[a])
    end
end

# Moves that were not precompiled: accumulate the pattern change of every touched cluster in the buffer's scratch space.
function _log_psi_diff_generic(x::AbstractConfig, move::AbstractMove, logψ::ClusterCPS{T}, Buffer::ClusterCPSBuffer) where {T}
    sites = affected_sites(move)
    dx = move_dx_before(move, x)
    (; pattern, delta, stamp, touched) = Buffer
    generation = Buffer.generation[] += 1
    ntouched = 0
    for (a, i) in enumerate(sites)
        δ = Int(dx[a])
        r = _site_range(logψ, i)
        (iszero(δ) || isempty(r)) && continue
        _digit(logψ, i, x[i] + δ)
        @inbounds for q in r
            c = logψ.site_cluster[q]
            if stamp[c] != generation
                stamp[c], delta[c] = generation, 0
                ntouched += 1
                touched[ntouched] = c
            end
            delta[c] += δ * logψ.site_stride[q]
        end
    end
    Δ = zero(T)
    @inbounds for n in 1:ntouched
        c = touched[n]
        old = logψ.cluster_offset[c] + pattern[c]
        Δ += logψ.params[old+delta[c]] - logψ.params[old]
    end
    return Δ
end

function post_move_affect!(Buffer::ClusterCPSBuffer, x::AbstractConfig, move::AbstractMove, logψ::ClusterCPS)
    j = _compiled_move_id(logψ, move)
    if !iszero(j)
        @inbounds for r in logψ.row_ptr[j]:logψ.row_ptr[j+1]-1
            c = logψ.row_cluster[r]
            Buffer.pattern[c] = _shifted(move, Buffer.pattern[c], logψ.row_shift[r])
        end
        return Buffer
    end
    sites = affected_sites(move)
    dx = move_dx_after(move, x)
    for (a, i) in enumerate(sites)
        δ = Int(dx[a])
        @inbounds for q in _site_range(logψ, i)
            Buffer.pattern[logψ.site_cluster[q]] += δ * logψ.site_stride[q]
        end
    end
    return Buffer
end
