# Helpers for ClusterCPS_tests.jl. Include after utils.jl.
using GreenFunctionMonteCarlo.LinearAlgebra

compareBuffers(B1::GFMC.ClusterCPSBuffer, B2::GFMC.ClusterCPSBuffer) = @test B1.pattern == B2.pattern

"Reference log ψ written directly from the definition, independent of the ClusterCPS data layout."
function reference_cps(x, clusters, tables, cluster_table; local_dim = 2, min_occupation = 0)
    s = 0.0
    for (c, cluster) in enumerate(clusters)
        I = sum((Int(x[site]) - min_occupation) * local_dim^(k - 1) for (k, site) in enumerate(cluster))
        s += tables[cluster_table[c]][I+1]
    end
    return s
end

"`n` random clusters on `sites`, with sizes drawn from `sizes` (sites within a cluster are distinct)."
random_clusters(rng, sites, n, sizes) = [shuffle(rng, sites)[1:rand(rng, sizes)] for _ in 1:n]

"Random tables, one per cluster size `L` in `sizes` and `copies` tables each; returns `(tables, cluster_table)`."
function random_tables(rng, clusters, T; local_dim = 2, copies = 2)
    sizes = sort(unique(length.(clusters)))
    tables = [randn(rng, T, local_dim^L) for L in sizes for _ in 1:copies]
    table_size = [L for L in sizes for _ in 1:copies]
    cluster_table = [rand(rng, findall(==(length(cl)), table_size)) for cl in clusters]
    return tables, cluster_table
end

"""
Random walk over applicable moves. At every step check, for every move of `H`, `log_psi_diff` against
`ψ(x′) - ψ(x)` (or -Inf if the move is not applicable), then apply a random applicable move, update the
buffer with `post_move_affect!` and compare it with a freshly computed one. Returns the number of checked moves.
"""
function check_random_walk(ψ, H, Hilbert, x0; nsteps = 100, rng = StableRNG(42), atol = 1e-10)
    x = copy(x0)
    B = GFMC.allocate_GWF_buffer(ψ, x)
    nchecked = 0
    for _ in 1:nsteps
        applicable = Int[]
        for (k, m) in enumerate(H.moves)
            Δ = GFMC.log_psi_diff(x, m, ψ, B, Hilbert)
            if GFMC.isapplicable(x, m, Hilbert)
                y = copy(x)
                apply!(y, m)
                @test Δ ≈ ψ(y) - ψ(x) atol = atol
                isfinite(Δ) && push!(applicable, k)
                nchecked += 1
            else
                @test Δ == -Inf
            end
        end
        isempty(applicable) && break
        m = H.moves[rand(rng, applicable)]
        apply!(x, m)
        GFMC.post_move_affect!(B, x, m, ψ)
        @test B.pattern == GFMC.allocate_GWF_buffer(ψ, x).pattern
    end
    return nchecked
end

"All configurations with occupations in `occupations` on `N` sites (optionally restricted by `filter`)."
product_basis(T, N, occupations; filter = x -> true) = Base.filter(filter, [collect(T, t) for t in Iterators.product(fill(occupations, N)...)])

"Dense matrix of `H` on `basis` (columns: initial states). Every move must map the basis into itself."
function dense_hamiltonian(H, Hilbert, basis)
    index = Dict(b => k for (k, b) in enumerate(basis))
    M = zeros(length(basis), length(basis))
    for (k, b) in enumerate(basis)
        x = BosonConfig(copy(b))
        M[k, k] += GFMC.get_diagonal(H)(x)
        for (m, w) in zip(H.moves, GFMC.get_offdiagonal_elements(H))
            GFMC.isapplicable(x, m, Hilbert) || continue
            y = copy(x)
            apply!(y, m)
            M[index[parent(y)], k] += w
        end
    end
    return M
end

"Ground state energy and (positive) ground state of `H` on `basis` by exact diagonalization."
function exact_ground_state(H, Hilbert, basis)
    M = dense_hamiltonian(H, Hilbert, basis)
    @assert M ≈ M'
    vals, vecs = eigen(Symmetric(M))
    ψ0 = vecs[:, 1]
    ψ0 .*= sign(sum(ψ0))
    @assert all(>(0), ψ0) "ground state not positive: Hamiltonian not ergodic on this basis?"
    return vals[1], ψ0
end

"Table of a single cluster covering `cluster` that equals `log ψ0` on `basis` (entries outside the basis stay `fill`)."
function exact_table(cluster, basis, ψ0; local_dim, min_occupation = 0, fill = 0.0)
    table = Base.fill(fill, local_dim^length(cluster))
    for (b, amplitude) in zip(basis, ψ0)
        I = sum((Int(b[site]) - min_occupation) * local_dim^(k - 1) for (k, site) in enumerate(cluster))
        table[I+1] = log(amplitude)
    end
    return table
end

"Sign-free soft-core boson model on a ring: constant-amplitude hopping -t, on-site energy U (n - n₀)(n - n₀ - 1)/2."
function soft_core_ring(N, n0, nmax; t = 1.0, U = 1.3)
    Hilbert = BosonHilbertSpace(N, OccupationNumberConstraint(n0, nmax))
    moves = Vector{Int8}[]
    for i in 1:N, j in (mod1(i - 1, N), mod1(i + 1, N))
        dx = zeros(Int8, N)
        dx[i], dx[j] = 1, -1
        push!(moves, dx)
    end
    Hxx = DiagOperator(x -> U * sum((Int(n) - n0) * (Int(n) - n0 - 1) / 2 for n in x))
    H = localOperator(moves, fill(-t, length(moves)), Hxx, Hilbert)
    return (; Hilbert, H)
end
