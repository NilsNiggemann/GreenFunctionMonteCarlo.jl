@testitem "ClusterCPS: pattern convention and evaluation" begin
    include("utils.jl")
    include("cps_utils.jl")

    @testset "hand-computed values" begin
        # pattern of cluster [3, 1] is x[3] + 2 x[1]: the first site is the least significant digit
        ψ = ClusterCPS([[3, 1]], [10.0, 20.0, 30.0, 40.0])
        @test ψ(Bool[0, 0, 0]) == 10.0
        @test ψ(Bool[0, 0, 1]) == 20.0
        @test ψ(Bool[1, 0, 0]) == 30.0
        @test ψ(Bool[1, 1, 1]) == 40.0

        # local_dim = 3, min_occupation = 1: pattern (x[2] - 1) + 3 (x[1] - 1)
        table = collect(1.0:9.0)
        ψ3 = ClusterCPS([[2, 1]], table; local_dim = 3, min_occupation = 1)
        @test ψ3(UInt8[1, 1]) == 1.0
        @test ψ3(UInt8[1, 3]) == 3.0
        @test ψ3(UInt8[3, 2]) == 8.0
        @test ψ3(UInt8[3, 3]) == 9.0

        # clusters of different sizes with different tables, plus one shared table
        ψmix = ClusterCPS([[1], [2, 3], [1]], [[1.0, 2.0], [0.0, 10.0, 20.0, 30.0]], [1, 2, 1])
        @test ψmix(Bool[1, 0, 1]) == 2.0 + 20.0 + 2.0
    end

    @testset "random clusters vs reference implementation" begin
        rng = StableRNG(7)
        for (local_dim, min_occupation, T) in ((2, 0, Bool), (3, 0, UInt8), (4, 2, UInt8))
            N = 9
            clusters = random_clusters(rng, 1:N-2, 12, 1:4)
            tables, cluster_table = random_tables(rng, clusters, Float64; local_dim)
            ψ = ClusterCPS(clusters, tables, cluster_table; local_dim, min_occupation)
            for _ in 1:50
                x = BosonConfig(T.(rand(rng, min_occupation:min_occupation+local_dim-1, N)))
                @test ψ(x) ≈ reference_cps(x, clusters, tables, cluster_table; local_dim, min_occupation) atol = 1e-12
                @test GFMC.allocate_GWF_buffer(ψ, x).pattern == [GFMC.cluster_pattern(x, ψ, c) for c in eachindex(clusters)]
            end
        end
    end

    @testset "equivalent constructions" begin
        rng = StableRNG(3)
        clusters = [[1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 1]]
        table = randn(rng, 8)
        shared = ClusterCPS(clusters, table)
        separate = ClusterCPS(clusters, [copy(table) for _ in clusters], 1:4)
        frommatrix = ClusterCPS(stack(clusters), table)
        # reversing the slot order of every cluster together with the table index gives the same state
        reversed_table = [table[1+sum(((I >> (3 - k)) & 1) << (k - 1) for k in 1:3)] for I in 0:7]
        reversed = ClusterCPS(reverse.(clusters), reversed_table)
        for _ in 1:20
            x = BosonConfig(rand(rng, Bool, 5))
            @test separate(x) ≈ shared(x)
            @test frommatrix(x) == shared(x)
            @test reversed(x) ≈ shared(x)
        end
        @test length(get_params(shared)) == 8
        @test length(get_params(separate)) == 32
    end

    @testset "parameters: shared tables, views, mutation" begin
        ψ = ClusterCPS([[1, 2], [2, 3], [1]], [zeros(4), zeros(2)], [1, 1, 2])
        x = BosonConfig(Bool[1, 1, 1])
        B = GFMC.allocate_GWF_buffer(ψ, x)
        tables = GFMC.get_tables(ψ)
        @test length(tables) == 2 && length.(tables) == [4, 2]
        tables[1][4] = 0.5          # pattern (1,1), used by both pair clusters
        tables[2][2] = 0.25
        @test get_params(ψ) == [0, 0, 0, 0.5, 0, 0.25]
        @test ψ(x) == 2 * 0.5 + 0.25
        # buffers stay valid after changing parameters: flipping site 2 removes both pair contributions
        Hilbert = BosonHilbertSpace(3, HardCoreConstraint())
        H = localOperator([Bool[0, 1, 0]], [-1.0], ZeroDiagOperator(), Hilbert)
        @test GFMC.log_psi_diff(x, H.moves[1], ψ, B, Hilbert) == -1.0
    end

    @testset "no clusters, uncovered sites, BitVector configs" begin
        Hilbert = BosonHilbertSpace(4, HardCoreConstraint())
        H = localOperator([Bool[1, 0, 0, 0], Bool[0, 0, 1, 1]], [-1.0, -1.0], ZeroDiagOperator(), Hilbert)
        x = BosonConfig(Bool[1, 0, 1, 0])

        empty = ClusterCPS(Vector{Int}[], Vector{Float64}[], Int[])
        @test empty(x) == 0.0
        B = GFMC.allocate_GWF_buffer(empty, x)
        @test all(m -> GFMC.log_psi_diff(x, m, empty, B, Hilbert) == 0.0, H.moves)

        # sites 3 and 4 are in no cluster; site 4 is beyond the largest cluster site
        ψ = ClusterCPS([[1, 2]], [0.0, 1.0, 2.0, 3.0])
        B = GFMC.allocate_GWF_buffer(ψ, x)
        @test GFMC.log_psi_diff(x, H.moves[2], ψ, B, Hilbert) == 0.0
        @test GFMC.log_psi_diff(x, H.moves[1], ψ, B, Hilbert) == -1.0
        @test check_random_walk(ψ, H, Hilbert, x; nsteps = 20) > 0

        xbits = BosonConfig(BitVector([1, 0, 1, 0]))
        @test ψ(xbits) == ψ(x)
    end

    @testset "Float32 tables" begin
        ψ = ClusterCPS([[1, 2], [2, 3]], Float32[0, 1, 2, 3])
        @test ψ isa ClusterCPS{Float32}
        Hilbert = BosonHilbertSpace(3, OccupationNumberConstraint(0, 1))
        H = localOperator([Int8[1, 0, 0], Int8[-1, 0, 0]], [-1.0, -1.0], ZeroDiagOperator(), Hilbert)
        x = BosonConfig(UInt8[1, 1, 0])
        B = GFMC.allocate_GWF_buffer(ψ, x)
        @test ψ(x) isa Float32
        @test GFMC.log_psi_diff(x, H.moves[2], ψ, B, Hilbert) === -1.0f0
        @test GFMC.log_psi_diff(x, H.moves[1], ψ, B, Hilbert) === -Inf32  # not applicable
    end

    @testset "integer tables are converted to floats" begin
        ψ = ClusterCPS([[1]], [0, 1])
        @test ψ isa ClusterCPS{Float64}
    end
end

@testitem "ClusterCPS: input validation" begin
    include("utils.jl")
    include("cps_utils.jl")

    @testset "constructor" begin
        @test_throws ArgumentError ClusterCPS([[1, 2]], zeros(3))                                   # wrong table length
        @test_throws ArgumentError ClusterCPS([[1, 2]], zeros(4); local_dim = 3)                    # needs 9
        @test_throws ArgumentError ClusterCPS([[1, 2], [1, 2, 3]], zeros(4))                        # shared table, sizes differ
        @test_throws ArgumentError ClusterCPS([[1, 1]], zeros(4))                                   # repeated site
        @test_throws ArgumentError ClusterCPS([Int[]], zeros(1))                                    # empty cluster
        @test_throws ArgumentError ClusterCPS([[0, 1]], zeros(4))                                   # site index < 1
        @test_throws ArgumentError ClusterCPS([[1]], zeros(1); local_dim = 1)
        @test_throws ArgumentError ClusterCPS([[1], [2]], [zeros(2)], [1, 2])                       # table index out of range
        @test_throws ArgumentError ClusterCPS([[1], [2]], [zeros(2)], [1])                          # cluster_table length
        @test_throws ArgumentError ClusterCPS([[1]], [0.0, NaN])
        @test_throws ArgumentError ClusterCPS([[1]], [0.0, Inf])
        @test ClusterCPS([[1]], [0.0, -Inf]) isa ClusterCPS                                         # -Inf (ψ = 0) is allowed
        @test ClusterCPS([[1], [2]], [zeros(2), zeros(4)], [1, 1]) isa ClusterCPS                   # unused tables are allowed
    end

    @testset "configurations outside the covered local states" begin
        ψ = ClusterCPS([[1, 2]], zeros(4))
        @test_throws ArgumentError ψ(UInt8[2, 0])
        @test_throws ArgumentError GFMC.allocate_GWF_buffer(ψ, BosonConfig(UInt8[0, 2]))
        @test_throws BoundsError GFMC.allocate_GWF_buffer(ψ, BosonConfig(UInt8[0]))                 # config shorter than clusters

        ψshift = ClusterCPS([[1]], zeros(2); min_occupation = 1)
        @test_throws ArgumentError ψshift(UInt8[0])

        # Hilbert space allows occupation 2, the CPS only 0:1: moves into 2 must fail loudly, not index a wrong entry
        Hilbert = BosonHilbertSpace(2, OccupationNumberConstraint(0, 2))
        H = localOperator([Int8[1, 0], Int8[0, 1]], [-1.0, -1.0], ZeroDiagOperator(), Hilbert)
        x = BosonConfig(UInt8[1, 0])
        B = GFMC.allocate_GWF_buffer(ψ, x)
        @test_throws ArgumentError GFMC.log_psi_diff(x, H.moves[1], ψ, B, Hilbert)
        @test GFMC.log_psi_diff(x, H.moves[2], ψ, B, Hilbert) == 0.0
        ψc = ClusterCPS([[1, 2]], zeros(4); moves = H.moves)
        @test GFMC._compiled_move_id(ψc, H.moves[1]) == 1
        @test_throws ArgumentError GFMC.log_psi_diff(x, H.moves[1], ψc, B, Hilbert)
        @test GFMC.log_psi_diff(x, H.moves[2], ψc, B, Hilbert) == 0.0
    end

    @testset "which Hilbert spaces make the target check redundant" begin
        covers(ψ, c) = GFMC._covers_hilbert_space(ψ, c)
        ψ2, ψ3, ψ3shift = ClusterCPS([[1]], zeros(2)), ClusterCPS([[1]], zeros(3); local_dim = 3), ClusterCPS([[1]], zeros(3); local_dim = 3, min_occupation = 1)
        @test covers(ψ2, HardCoreConstraint()) && covers(ψ3, HardCoreConstraint()) && !covers(ψ3shift, HardCoreConstraint())
        @test covers(ψ2, OccupationNumberConstraint(0, 1)) && !covers(ψ2, OccupationNumberConstraint(0, 2))
        @test covers(ψ3, OccupationNumberConstraint(0, 2)) && covers(ψ3, OccupationNumberConstraint(1, 2))
        @test covers(ψ3shift, OccupationNumberConstraint(1, 3)) && !covers(ψ3shift, OccupationNumberConstraint(0, 2)) && !covers(ψ3shift, OccupationNumberConstraint(1, 4))
        @test !covers(ψ2, GFMC.NoConstraint())
    end

    @testset "moves that cannot be precompiled" begin
        Hflip = localOperator([Bool[1, 0]], [-1.0], ZeroDiagOperator(), BosonHilbertSpace(2, HardCoreConstraint()))
        @test_throws ArgumentError ClusterCPS([[1, 2]], zeros(9); local_dim = 3, moves = Hflip.moves)
        @test_throws ArgumentError ClusterCPS([[1, 2]], zeros(4); min_occupation = 1, moves = Hflip.moves)
        struct StateDependentMove <: GFMC.AbstractMove
            inds::Vector{Int}
        end
        GFMC.affected_sites(m::StateDependentMove) = m.inds
        @test_throws ArgumentError ClusterCPS([[1, 2]], zeros(4); moves = [StateDependentMove([1])])
    end
end

@testitem "ClusterCPS: log_psi_diff and buffer updates along random walks" begin
    include("utils.jl")
    include("cps_utils.jl")
    rng = StableRNG(11)
    N = 10
    covered = 1:N-2                       # sites N-1, N are in no cluster, but moves touch them

    @testset "hard-core, flip moves of $(isnothing(len) ? "variable" : "fixed") length, $T, compiled = $compiled" for len in (3, nothing), T in (Float64, Float32), compiled in (false, true)
        clusters = random_clusters(rng, covered, 14, 1:5)
        tables, cluster_table = random_tables(rng, clusters, T)
        (Hilbert, H) = getExample(N, 30, rng, len, HardCoreConstraint())
        ψ = ClusterCPS(clusters, tables, cluster_table; moves = compiled ? H.moves : nothing)
        x = BosonConfig(rand(rng, Bool, N))
        @test check_random_walk(ψ, H, Hilbert, x; atol = T == Float32 ? 1e-4 : 1e-10) > 500
        TestWFRatio(ψ, x, H, Hilbert; tol = T == Float32 ? 1e-4 : 1e-10)
        testPostMove(ψ, x, H, Hilbert)
    end

    @testset "soft-core, local_dim = $local_dim, min_occupation = $n0, compiled = $compiled" for (local_dim, n0) in ((3, 0), (3, 1), (4, 0)), compiled in (false, true)
        nmax = n0 + local_dim - 1
        Hilbert = BosonHilbertSpace(N, OccupationNumberConstraint(n0, nmax))
        moves = Vector{Int8}[]
        for _ in 1:40      # hopping by one or two particles, and creation/annihilation of up to local_dim - 1
            dx = zeros(Int8, N)
            i, j = shuffle(rng, 1:N)[1:2]
            if rand(rng, Bool)
                k = rand(rng, 1:2)
                dx[i], dx[j] = k, -k
            else
                dx[i] = rand(rng, [-(local_dim - 1):-1; 1:local_dim-1])
            end
            push!(moves, dx)
        end
        H = localOperator(moves, -rand(rng, length(moves)), ZeroDiagOperator(), Hilbert)
        clusters = random_clusters(rng, covered, 10, 1:3)
        tables, cluster_table = random_tables(rng, clusters, Float64; local_dim)
        ψ = ClusterCPS(clusters, tables, cluster_table; local_dim, min_occupation = n0, moves = compiled ? H.moves : nothing)
        x = BosonConfig(UInt8.(rand(rng, n0:nmax, N)))
        @test check_random_walk(ψ, H, Hilbert, x) > 500
    end

    @testset "move inside a single cluster, and across many overlapping clusters, compiled = $compiled" for compiled in (false, true)
        Hilbert = BosonHilbertSpace(6, HardCoreConstraint())
        # all affected sites in one cluster (pattern changes by several digits at once)
        H = localOperator([Bool[1, 1, 1, 0, 0, 0], Bool[0, 1, 0, 1, 0, 1], Bool[1, 1, 1, 1, 1, 1]], -ones(3), ZeroDiagOperator(), Hilbert)
        clusters = [[1, 2, 3], [3, 2, 1], [2, 4, 6], [1, 2, 3, 4, 5, 6], [6], [5, 4]]   # identical site sets in different slot orders
        tables, cluster_table = random_tables(rng, clusters, Float64)
        ψ = ClusterCPS(clusters, tables, cluster_table; moves = compiled ? H.moves : nothing)
        x = BosonConfig(rand(rng, Bool, 6))
        @test check_random_walk(ψ, H, Hilbert, x) > 0
    end

    @testset "scratch space does not leak between calls" begin
        clusters = random_clusters(rng, 1:N, 12, 2:4)
        tables, cluster_table = random_tables(rng, clusters, Float64)
        ψ = ClusterCPS(clusters, tables, cluster_table)
        (Hilbert, H) = getExample(N, 25, rng, nothing, HardCoreConstraint())
        x = BosonConfig(rand(rng, Bool, N))
        B = GFMC.allocate_GWF_buffer(ψ, x)
        forward = [GFMC.log_psi_diff(x, m, ψ, B, Hilbert) for m in H.moves]
        backward = reverse([GFMC.log_psi_diff(x, m, ψ, B, Hilbert) for m in reverse(H.moves)])
        repeated = [(GFMC.log_psi_diff(x, m, ψ, B, Hilbert); GFMC.log_psi_diff(x, m, ψ, B, Hilbert)) for m in H.moves]
        @test forward == backward == repeated
        fresh = [GFMC.log_psi_diff(x, m, ψ, GFMC.allocate_GWF_buffer(ψ, x), Hilbert) for m in H.moves]
        @test forward == fresh
    end

    @testset "setBuffer! copies the walker state" begin
        clusters = random_clusters(rng, 1:N, 8, 1:3)
        tables, cluster_table = random_tables(rng, clusters, Float64)
        ψ = ClusterCPS(clusters, tables, cluster_table)
        x, y = BosonConfig(rand(rng, Bool, N)), BosonConfig(rand(rng, Bool, N))
        A, B = GFMC.allocate_GWF_buffer(ψ, x), GFMC.allocate_GWF_buffer(ψ, y)
        GFMC.log_psi_diff(x, GFMC.FlipMove(SC.SmallVector{2,Int}((1, 2))), ψ, A, BosonHilbertSpace(N, HardCoreConstraint()))
        GFMC.setBuffer!(A, B)
        @test A.pattern == B.pattern && A.pattern !== B.pattern
        (Hilbert, H) = getExample(N, 10, rng, 2, HardCoreConstraint())
        @test [GFMC.log_psi_diff(y, m, ψ, A, Hilbert) for m in H.moves] == [GFMC.log_psi_diff(y, m, ψ, B, Hilbert) for m in H.moves]
    end
end

@testitem "ClusterCPS: -Inf entries forbid patterns" begin
    include("utils.jl")
    include("cps_utils.jl")
    # nearest-neighbour exclusion on an open chain: pattern (1,1) of every bond has ψ = 0
    N, NWalkers, NSteps = 8, 20, 50
    Hilbert = BosonHilbertSpace(N, HardCoreConstraint())
    H = localOperator(collect(eachcol(Bool.(I(N)))), -ones(N), DiagOperator(x -> 0.3 * sum(x)), Hilbert)
    @testset "compiled = $(!isnothing(moves))" for moves in (nothing, H.moves)
    ψ = ClusterCPS([[i, i + 1] for i in 1:N-1], [0.0, 0.1, 0.1, -Inf]; moves)
    x = BosonConfig(Bool[1, 0, 0, 1, 0, 1, 0, 0])
    B = GFMC.allocate_GWF_buffer(ψ, x)
    @test GFMC.log_psi_diff(x, H.moves[2], ψ, B, Hilbert) == -Inf     # would occupy sites 1 and 2
    @test isfinite(GFMC.log_psi_diff(x, H.moves[8], ψ, B, Hilbert))
    @test isfinite(GFMC.getLocalEnergy(x, H, ψ, Hilbert))
    @test check_random_walk(ψ, H, Hilbert, x) > 0

    prob = GFMCProblem(x, NWalkers, ContinuousTimePropagator(0.5); logψ = ψ, H, Hilbert, parallelization = GFMC.SingleThreaded())
    Obs = ConfigObserver(x, NSteps, NWalkers)
    runGFMC!(prob, Obs, NSteps; rng = StableRNG(5), logger = nothing)
    configs = GFMC.getObs(Obs).SaveConfigs
    @test !any(configs[i, α, n] && configs[i+1, α, n] for i in 1:N-1, α in 1:NWalkers, n in 1:NSteps)
    @test all(isfinite, GFMC.getObs(Obs).energies)
    end
end

@testitem "ClusterCPS: pair and single-site clusters reproduce a Jastrow factor" begin
    include("utils.jl")
    include("cps_utils.jl")
    NSites, NWalkers, NSteps = 7, 20, 30
    rng = StableRNG(21)
    m = randn(rng, NSites)
    v = randn(rng, NSites, NSites)
    v = v + v'
    jastrow = Jastrow(m, v)

    # log ψ = Σ_i m_i n_i + ½ Σ_ij v_ij n_i n_j, with n_i² = n_i for hard-core bosons
    singles = [[i] for i in 1:NSites]
    pairs = [[i, j] for i in 1:NSites for j in i+1:NSites]
    tables = [[[0.0, m[i] + v[i, i] / 2] for i in 1:NSites]; [[0.0, 0.0, 0.0, v[i, j]] for (i, j) in pairs]]
    ψ = ClusterCPS([singles; pairs], tables, eachindex(tables))

    (Hilbert, H) = getExample(NSites, 12, rng, nothing, HardCoreConstraint())
    x = BosonConfig(rand(rng, Bool, NSites))
    Bj, Bc = GFMC.allocate_GWF_buffer(jastrow, x), GFMC.allocate_GWF_buffer(ψ, x)
    @test ψ(x) ≈ jastrow(x)
    @test [GFMC.log_psi_diff(x, mv, ψ, Bc, Hilbert) for mv in H.moves] ≈ [GFMC.log_psi_diff(x, mv, jastrow, Bj, Hilbert) for mv in H.moves]

    run_gfmc(logψ) = begin
        prob = GFMCProblem(x, NWalkers, ContinuousTimePropagator(1.0); logψ, H, Hilbert, parallelization = GFMC.SingleThreaded())
        Obs = ConfigObserver(x, NSteps, NWalkers)
        runGFMC!(prob, Obs, NSteps; rng = StableRNG(1234), logger = nothing)
        GFMC.getObs(Obs)
    end
    obs_jastrow = run_gfmc(jastrow)
    for logψ in (ψ, ClusterCPS([singles; pairs], tables, eachindex(tables); moves = H.moves))
        obs_cps = run_gfmc(logψ)
        @test obs_cps.SaveConfigs == obs_jastrow.SaveConfigs
        @test obs_cps.energies ≈ obs_jastrow.energies atol = 1e-10
    end
end

@testitem "ClusterCPS: GFMC runs match the unbuffered NaiveFunction" begin
    include("utils.jl")
    include("cps_utils.jl")
    NSites = 9
    rng = StableRNG(99)
    clusters = random_clusters(rng, 1:NSites, 9, 2:4)
    tables, cluster_table = random_tables(rng, clusters, Float64)
    tables .*= 0.3
    ψ = ClusterCPS(clusters, tables, cluster_table)

    @testset "hard-core, $(isnothing(len) ? "variable" : "fixed") move length" for len in (2, nothing)
        (Hilbert, H) = getExample(NSites, 15, rng, len, HardCoreConstraint())
        x = BosonConfig(rand(rng, Bool, NSites))
        testRun(x, ψ, H, Hilbert; NWalkers = 20, NSteps = 20)
        testRun(x, ClusterCPS(clusters, tables, cluster_table; moves = H.moves), H, Hilbert; NWalkers = 20, NSteps = 20)
    end

    @testset "soft-core ring" begin
        (; Hilbert, H) = soft_core_ring(NSites, 0, 2)
        tables3, cluster_table3 = random_tables(rng, clusters, Float64; local_dim = 3)
        x = BosonConfig(UInt8.(rand(rng, 0:2, NSites)))
        for moves in (nothing, H.moves)
            ψ3 = ClusterCPS(clusters, 0.3 .* tables3, cluster_table3; local_dim = 3, moves)
            testRun(x, ψ3, H, Hilbert; NWalkers = 20, NSteps = 20)
        end
    end

    @testset "multithreaded, $(nameof(typeof(par)))" for par in (GFMC.MultiThreaded(4), GFMC.BatchMultiThreaded(4))
        (Hilbert, H) = getExample(NSites, 15, rng, nothing, HardCoreConstraint())
        x = BosonConfig(rand(rng, Bool, NSites))
        results = map((NaiveFunction(ψ), ψ, ClusterCPS(clusters, tables, cluster_table; moves = H.moves))) do logψ
            prob = GFMCProblem(x, 32, ContinuousTimePropagator(1.0); logψ, H, Hilbert, parallelization = par)
            Obs = ConfigObserver(x, 20, 32)
            runGFMC!(prob, Obs, 20; rng = StableRNG(8), logger = nothing)
            GFMC.getObs(Obs)
        end
        @test results[1].SaveConfigs == results[2].SaveConfigs == results[3].SaveConfigs
        @test results[1].energies ≈ results[2].energies atol = 1e-10
        @test results[1].energies ≈ results[3].energies atol = 1e-10
    end
end

@testitem "ClusterCPS: exact ground state as a single cluster gives zero-variance energies" begin
    include("utils.jl")
    include("cps_utils.jl")

    # For the generic and the precompiled guide: E_loc(x) = E0 for every x in `basis`, and every GFMC step returns E0
    function test_zero_variance(clusters, tables, cluster_table, H, Hilbert, basis, E0, x0; NWalkers = 10, NSteps = 30, kwargs...)
        for moves in (nothing, H.moves)
            ψ = ClusterCPS(clusters, tables, cluster_table; moves, kwargs...)
            @test all(b -> isapprox(GFMC.getLocalEnergy(BosonConfig(copy(b)), H, ψ, Hilbert), E0; atol = 1e-9), basis)
            prob = GFMCProblem(x0, NWalkers, ContinuousTimePropagator(0.3); logψ = ψ, H, Hilbert)
            Obs = BasicObserver(NSteps, NWalkers)
            runGFMC!(prob, Obs, NSteps; rng = StableRNG(3), logger = nothing)
            @test all(E -> isapprox(E, E0; atol = 1e-9), Obs.energies)
        end
    end

    @testset "transverse-field Ising chain, sites in permuted slot order" begin
        N = 6
        (; Hilbert, H) = getMinimalExample(N, 0.7, 1.0)
        basis = product_basis(Bool, N, 0:1)
        E0, ψ0 = exact_ground_state(H, Hilbert, basis)
        cluster = shuffle(StableRNG(1), 1:N)
        test_zero_variance([cluster], [exact_table(cluster, basis, ψ0; local_dim = 2)], [1], H, Hilbert, basis, E0, BosonConfig(Bool[1, 0, 1, 1, 0, 0]))
    end

    @testset "two decoupled chains: product of two cluster tables" begin
        # sites 1:3 and 4:6 form two independent TFI chains, so ψ0 = ψa ⊗ ψb exactly
        Na = 3
        Hilbert = BosonHilbertSpace(2Na, HardCoreConstraint())
        Hxx = DiagOperator(x -> sum(σz(i, x) * σz(i + 1, x) for i in (1, 2, 4, 5)))
        H = localOperator(collect(eachcol(Bool.(I(2Na)))), -0.8 .* ones(2Na), Hxx, Hilbert)
        Hilbert_a = BosonHilbertSpace(Na, HardCoreConstraint())
        Ha = localOperator(collect(eachcol(Bool.(I(Na)))), -0.8 .* ones(Na), DiagOperator(x -> sum(σz(i, x) * σz(i + 1, x) for i in 1:Na-1)), Hilbert_a)
        basis_a = product_basis(Bool, Na, 0:1)
        Ea, ψa = exact_ground_state(Ha, Hilbert_a, basis_a)
        table = exact_table(1:Na, basis_a, ψa; local_dim = 2)
        test_zero_variance([[1, 2, 3], [4, 5, 6]], [table, copy(table)], [1, 2], H, Hilbert, product_basis(Bool, 2Na, 0:1), 2Ea, BosonConfig(Bool[0, 1, 0, 1, 1, 0]))
    end

    @testset "soft-core ring, min_occupation = $n0" for n0 in (0, 1)
        N, nmax, Npart = 4, n0 + 2, 4
        (; Hilbert, H) = soft_core_ring(N, n0, nmax)
        basis = product_basis(UInt8, N, n0:nmax; filter = b -> sum(Int, b) - N * n0 == Npart)
        E0, ψ0 = exact_ground_state(H, Hilbert, basis)
        cluster = [2, 4, 1, 3]
        table = exact_table(cluster, basis, ψ0; local_dim = 3, min_occupation = n0, fill = -Inf)
        test_zero_variance([cluster], [table], [1], H, Hilbert, basis, E0, BosonConfig(UInt8.(n0 .+ [1, 1, 1, 1])); local_dim = 3, min_occupation = n0)
    end
end

@testitem "ClusterCPS: precompiled moves" begin
    include("utils.jl")
    include("cps_utils.jl")
    rng = StableRNG(5)
    N = 8
    clusters = random_clusters(rng, 1:N-1, 10, 1:4)
    tables, cluster_table = random_tables(rng, clusters, Float64)

    @testset "every move of H is found, and agrees with the generic path" begin
        for (Hilbert, H) in (getExample(N, 20, rng, 3, HardCoreConstraint()), soft_core_ring(N, 0, 1))
            generic = ClusterCPS(clusters, tables, cluster_table)
            compiled = ClusterCPS(clusters, tables, cluster_table; moves = H.moves)
            @test all(m -> GFMC._compiled_move_id(compiled, m) > 0, H.moves)
            @test all(m -> GFMC._compiled_move_id(generic, m) == 0, H.moves)
            for _ in 1:20
                x = BosonConfig(Hilbert)
                x .= rand(rng, Bool, N)
                Bg, Bc = GFMC.allocate_GWF_buffer(generic, x), GFMC.allocate_GWF_buffer(compiled, x)
                @test [GFMC.log_psi_diff(x, m, compiled, Bc, Hilbert) for m in H.moves] ≈ [GFMC.log_psi_diff(x, m, generic, Bg, Hilbert) for m in H.moves] atol = 1e-12
            end
        end
    end

    @testset "lookup is by identity of the move, foreign moves use the generic path" begin
        Hilbert = BosonHilbertSpace(N, OccupationNumberConstraint(0, 1))
        function hop(i, j)
            dx = zeros(Int8, N)
            dx[i], dx[j] = 1, -1
            return dx
        end
        H = localOperator([hop(1, 2), hop(3, 4)], [-1.0, -1.0], ZeroDiagOperator(), Hilbert)
        ψ = ClusterCPS(clusters, tables, cluster_table; moves = H.moves)
        generic = ClusterCPS(clusters, tables, cluster_table)
        H_other = localOperator([hop(2, 1), hop(1, 3), hop(3, 4)], -ones(3), ZeroDiagOperator(), Hilbert)
        @test GFMC._compiled_move_id(ψ, H_other.moves[1]) == 0       # same sites, other values
        @test GFMC._compiled_move_id(ψ, H_other.moves[2]) == 0       # same first site, other sites
        @test GFMC._compiled_move_id(ψ, H_other.moves[3]) == 2       # equal (isbits) move from another operator
        flip = GFMC.FlipMove(H.moves[1].inds)
        @test GFMC._compiled_move_id(ψ, flip) == 0                  # same sites, other move type
        x = BosonConfig(UInt8[0, 1, 1, 0, 1, 0, 1, 0])
        B, Bg = GFMC.allocate_GWF_buffer(ψ, x), GFMC.allocate_GWF_buffer(generic, x)
        for m in H_other.moves
            @test GFMC.log_psi_diff(x, m, ψ, B, Hilbert) ≈ GFMC.log_psi_diff(x, m, generic, Bg, Hilbert) atol = 1e-12
        end
    end

    @testset "duplicate, empty and heterogeneous move lists" begin
        Hilbert = BosonHilbertSpace(N, HardCoreConstraint())
        H = localOperator([Bool[1, 1, 0, 0, 0, 0, 0, 0], zeros(Bool, N), Bool[1, 1, 0, 0, 0, 0, 0, 0], Bool[0, 0, 0, 0, 0, 0, 0, 1]], -ones(4), ZeroDiagOperator(), Hilbert)
        ψ = ClusterCPS(clusters, tables, cluster_table; moves = H.moves)   # site 8 is in no cluster
        @test GFMC._compiled_move_id(ψ, H.moves[2]) == 0                    # no affected sites
        @test GFMC._compiled_move_id(ψ, H.moves[3]) == 1                    # first of the duplicates
        @test check_random_walk(ψ, H, Hilbert, BosonConfig(rand(rng, Bool, N))) > 0

        mixed = GFMC.AbstractMove[GFMC.FlipMove(SC.SmallVector{2,Int}((1, 2))), GFMC.SparseMove([3, 4], Int8[1, -1])]
        ψmixed = ClusterCPS(clusters, tables, cluster_table; moves = mixed)
        @test GFMC._compiled_move_id(ψmixed, mixed[1]) == 1
        @test GFMC._compiled_move_id(ψmixed, mixed[2]) == 2               # non-isbits move: same object
        x = BosonConfig(UInt8[1, 0, 0, 1, 0, 0, 1, 0])
        Bm, Bg = GFMC.allocate_GWF_buffer(ψmixed, x), GFMC.allocate_GWF_buffer(ClusterCPS(clusters, tables, cluster_table), x)
        H_occ = BosonHilbertSpace(N, OccupationNumberConstraint(0, 1))
        @test GFMC.log_psi_diff(x, mixed[2], ψmixed, Bm, H_occ) ≈ GFMC.log_psi_diff(x, mixed[2], ClusterCPS(clusters, tables, cluster_table), Bg, H_occ) atol = 1e-12
    end
end
