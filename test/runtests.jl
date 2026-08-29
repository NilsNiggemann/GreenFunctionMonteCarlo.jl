using TestItemRunner, Test

@testitem "Bosonic Configuration Tests" begin
    include("utils.jl")
    Hilbert = BosonHilbertSpace(10, OccupationNumberConstraint(0, 1))
    RNG = StableRNG(1234)
    config = BosonConfig(Hilbert)
    
    @testset "Int8 Variables" begin
        @test size(config) === (10,)
        @test config == BosonConfig(zeros(Int8,10))
        @test fulfills_constraint(config, Hilbert)
        rand!(RNG,config)
        @test !fulfills_constraint(config, Hilbert)
    end

    config1 = BosonConfig(zeros(Bool,10))
    config2 = BosonConfig(BitVector(zeros(Bool,10)))
    
    for config in (config1, config2)
        @testset "Bool Variables" begin
            @test size(config) === (10,)
            @test fulfills_constraint(config, Hilbert)
            rand!(RNG,config)
            @test fulfills_constraint(config, Hilbert)
        end
    end

    @testset "HardCoreConstraint" begin
        Hilbert = BosonHilbertSpace(10, HardCoreConstraint())
        config = BosonConfig(Hilbert)
        @test fulfills_constraint(config, Hilbert)
        rand!(RNG,config)
        @test fulfills_constraint(config, Hilbert)
    end

end

@testitem "AbstractConfigs" begin
    include("utils.jl")
    struct TestConfig{T} <: AbstractConfig{T,1}
        data::Matrix{T}
    end
    Base.parent(x::TestConfig) = x.data

    x = TestConfig(zeros(Int,10,2))
    @test getindex(x,1) == 0
    setindex!(x,1,1)
    @test x[1] == 1
    @test getindex(x,1,1) == 1

    @test sum(x) == 1

end
@testitem "Sparse Move Tests" begin
    include("utils.jl")
    moves = [
        Int8[0,-1,1],
        Int8[0,1,-1],
        Int8[1,0,1],
        Int8[-1,0,-1],
    ]

    moves_binary = [
        Bool[0,1,1],
        Bool[1,0,1]
    ]

    weights = -[0.5, 0.5, 0.2,0.2]
    weights_bin = -[0.5, 0.2]
    
    Hilbert = BosonHilbertSpace(3, HardCoreConstraint())

    diag = ZeroDiagOperator()
    local_op = localOperator(moves, weights, diag, Hilbert)

    config = BosonConfig(Int8[1,1,0])
    config_Bin = BosonConfig(Bool[1,1,0])

    @testset "SparseMove" begin
        move = GFMC.get_move(local_op, 1)
        GFMC.apply!(config, move)
        @test config == [1,0,1]
    end

    @testset "SparseMoveBinary" begin
        local_op_binary = localOperator(moves_binary, weights_bin, diag, Hilbert)
        move_binary = GFMC.get_move(local_op_binary, 1)
        GFMC.apply!(config_Bin, move_binary)
        @test config_Bin == Bool[1,0,1]
    end
end

##
@testitem "Diagonal Operator Tests" begin
    include("utils.jl")
    Nsites = 5
    Hilbert = BosonHilbertSpace(Nsites, HardCoreConstraint())
    x = BosonConfig(Hilbert)
    
    RNG = StableRNG(1234)
    rand!(RNG,x)

    d = rand(RNG, Nsites)
    v = rand(RNG, Nsites, Nsites)

    H1 = OneBodyDiagOperator(d)
    H2 = TwoBodyDiagOperator(v)

    Hsum = H1 + H2

    @testset "Diagonal Operators" begin
        @test Hsum isa DiagOperatorSum
        @test H1(x) == d'*x
        @test H2(x) == x'v*x
        @test Hsum(x) == H1(x)+H2(x)
    end
end
##
@testitem "Bosonic Walker Ensemble Tests" begin
    include("utils.jl")
    Hilbert = BosonHilbertSpace(3, HardCoreConstraint())
    config = BosonConfig(Hilbert)
    logψ = GFMC.EqualWeightSuperposition()
    
    CT = ContinuousTimePropagator(1.)
    
    H = localOperator(
        [
            Bool[0,1,1],
            Bool[0,0,1],
            Bool[0,1,0],
            Bool[1,1,0],
            ]
    , -[0.5, 0.3, 0.2 ,0.4], ZeroDiagOperator(), Hilbert)
    
    NumWalkers = 8
    ensemble = GFMC.allocate_walkerEnsemble(config,logψ,NumWalkers,H)
    @testset "FlipMoves" begin
        move = H.moves[1]
        @test move === GFMC.FlipMove(SC.SmallVector{2,Int}((2,3)))

        GFMC.move_dx_before(H.moves[1],config) === SC.SmallVector{2,Bool}((1,1))
        @test iszero(config)
        
        GFMC.apply!(config, move)
        @test config == Bool[0,1,1]
        GFMC.move_dx_after(H.moves[1],config) === SC.SmallVector{2,Bool}((1,1))

        GFMC.apply_inverse!(config, move)
        

    end
    @testset "Walker Ensemble Construction" begin
        @test length(ensemble.Configs) == NumWalkers
        @test length(ensemble.WalkerWeights) == NumWalkers
        @test length(ensemble.MoveWeights) == NumWalkers
        @test length(ensemble.Buffers) == NumWalkers
        @test all(length.(ensemble.MoveWeights) .== 4)
    end
    
    move_weights = GFMC.getMoveWeights(ensemble,1)
    
    GFMC.get_markov_weights!(move_weights, config, H, logψ, Hilbert, GFMC.getBuffer(ensemble, 1))
    operator_weights = GFMC.get_offdiagonal_elements(H)
    
    @testset "Markov Weights" begin
        @test all(>=(0), move_weights)
        @test move_weights == -operator_weights
    end
    
    RNG = [StableRNG(1234)]
    GFMC.propagateWalkers!(ensemble, H, logψ, Hilbert, CT, GFMC.SingleThreaded(),RNG)
    
    AllConfs = stack(ensemble.Configs)
    @testset "Continuous Time Propagation" begin
        @test AllConfs == Bool[0 0 0 0 1 0 0 0; 1 0 1 0 1 1 0 1; 0 1 1 0 0 0 0 1]
        GFMC.propagateWalkers!(ensemble, H, logψ, Hilbert, CT, GFMC.SingleThreaded(),RNG)

        @test stack(ensemble.Configs) != AllConfs
    end

end
##
@testitem "Discrete Propagator Walker Ensemble Tests" begin
    include("utils.jl")
    Hilbert = BosonHilbertSpace(3, HardCoreConstraint())
    config = BosonConfig(Hilbert)
    logψ = GFMC.EqualWeightSuperposition()

    DT = DiscretePropagator(1., 3)

    H = localOperator(
        [
            Bool[0,1,1],
            Bool[0,0,1],
            Bool[0,1,0],
            Bool[1,1,0],
            ]
    , -[0.5, 0.3, 0.2 ,0.4], ZeroDiagOperator(), Hilbert)

    NumWalkers = 8
    ensemble = GFMC.allocate_walkerEnsemble(config,logψ,NumWalkers,H)

    RNG = [StableRNG(1234)]
    GFMC.propagateWalkers!(ensemble, H, logψ, Hilbert, DT, GFMC.SingleThreaded(),RNG)

    @testset "Weights and Energies" begin
        @test all(isfinite, GFMC.getWalkerWeights(ensemble))
        @test all(>=(0), GFMC.getWalkerWeights(ensemble))
        @test all(isfinite, GFMC.getLocalEnergies(ensemble))
    end

    AllConfs = stack(ensemble.Configs)
    @testset "Discrete Time Propagation" begin
        @test AllConfs != zeros(Bool,3,NumWalkers)
    end

    @testset "Reproducibility" begin
        ensemble2 = GFMC.allocate_walkerEnsemble(config,logψ,NumWalkers,H)
        RNG2 = [StableRNG(1234)]
        GFMC.propagateWalkers!(ensemble2, H, logψ, Hilbert, DT, GFMC.SingleThreaded(),RNG2)

        @test stack(ensemble2.Configs) == AllConfs
        @test GFMC.getWalkerWeights(ensemble2) == GFMC.getWalkerWeights(ensemble)
        @test GFMC.getLocalEnergies(ensemble2) == GFMC.getLocalEnergies(ensemble)
    end
end
##


@testitem "main Usage" begin
    include("utils.jl")
    
    (;Hilbert,H) = getExampleHardcore(3,4,StableRNG(1234))

    config = BosonConfig(Hilbert)
    RNG = StableRNG(1234)
    rand!(RNG,config)
    logψ = GFMC.EqualWeightSuperposition()
    NWalkers = 10
    CT = ContinuousTimePropagator(1.)
    

    prob = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert)

    @testset "runGFMC" begin
        obs = NoObserver()
        runGFMC!(prob, obs, 10; rng= RNG)
        AllConfs = stack(prob.Walkers.Configs)
        @test AllConfs != zeros(Bool,3,NWalkers)
    end

    outfile = tempname()

    NSteps = 5

    ConfSaverFile = ConfigObserver(outfile, config, NSteps, NWalkers)

    @testset "ConfigObserver" begin

        @test isfile(outfile)
        runGFMC!(prob, ConfSaverFile, NSteps; rng= RNG)

        GFMC.HDF5.h5open(outfile, "r") do file
            @test haskey(file, "SaveConfigs")

            SaveConfigs = read(file["SaveConfigs"])

            @test haskey(file, "TotalWeights")
            TotalWeights = read(file["TotalWeights"])
            @test haskey(file, "energies")
            energies = read(file["energies"])
            @test haskey(file, "reconfigurationTable")
            reconfigurationTable = read(file["reconfigurationTable"])

            testSaveConf(SaveConfigs,TotalWeights,energies,reconfigurationTable,3,NWalkers,NSteps)
        end
    end
end

@testitem "Accumulator TFI" begin
    include("utils.jl")
    using GreenFunctionMonteCarlo.Statistics:mean

    function E_critPoint_exact(L)
        return 1 - csc(pi / (2 * (2 * L + 1)))
    end

    using GreenFunctionMonteCarlo.LinearAlgebra
    NSites = 2
    NSteps = 1500
    mProj = 20
    NWalkers = 2

    RNG = StableRNG(1234)

    Hilbert = BosonHilbertSpace(NSites, HardCoreConstraint())
    moves = eachcol(Bool.(I(NSites))) # each move flips a single spin
    offdiagElements = -ones(NSites)
    Hxx = DiagOperator(x-> sum(σz(i, x) * σz(i + 1, x) for i in eachindex(x)[1:end-1]))

    H = localOperator(moves, offdiagElements, Hxx, Hilbert)

    config = BosonConfig(Hilbert)
    rand!(RNG,config)
    logψ = GFMC.EqualWeightSuperposition()
    CT = ContinuousTimePropagator(0.1)
    
    prob = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert)

    # Estimate weights for the continuous time propagator

    weight_normalization, w_avg_estimate = GFMC.estimate_weights_continuousTime!(prob;Nepochs=3,Nsamples=30,mProj = 20,rng = RNG)
    CT = ContinuousTimePropagator(0.1, w_avg_estimate)
    
    outfile = tempname()

    BObs = GFMC.BasicObserver(outfile, NSteps, NWalkers)
    CObs = GFMC.ConfigurationObserver(outfile, config, NSteps, NWalkers)
    
    BasicAccumulatorFile = GFMC.BasicAccumulator(outfile, mProj, NWalkers; weight_normalization)
    ObsAccumulatorFile = GFMC.ObservableAccumulator(outfile,OccupationNumber(NSites),BasicAccumulatorFile, mProj, NWalkers, Threads.nthreads())

    Observer = GFMC.CombinedObserver((BObs, CObs,BasicAccumulatorFile, ObsAccumulatorFile))

    runGFMC!(prob, Observer, NSteps; rng = RNG)

    Energy = GFMC.getEnergies(BObs.TotalWeights, BObs.energies, mProj)
    # Energy_direct = BasicAccumulatorFile.en_numerator ./ BasicAccumulatorFile.Gnp_denominator .*NSites
    Energy_direct = mean(GFMC.get_energy_from_accumulator_bunching(BasicAccumulatorFile, 1)) .* NSites
    @testset "BasicAccumulator" begin

        @test isfile(outfile)

        GFMC.HDF5.h5open(outfile, "r") do file
            @test haskey(file, "Gnp_denominator")

            Gnp_denominator = read(file["Gnp_denominator"])

            @test haskey(file, "en_numerator")
            en_numerator = read(file["en_numerator"])
            @test !iszero(Gnp_denominator)
            @test !iszero(en_numerator)

            @test Energy ≈ Energy_direct atol = 1e-10
            @test Energy[end÷2] ≈ E_critPoint_exact(NSites) rtol = 3e-2
        end
    end

    Gnps = GFMC.precomputeNormalizedAccWeight(BObs.TotalWeights,2mProj)

    n_avg = stack(getObs_diagonal(Gnps,CObs.SaveConfigs,BObs.reconfigurationTable,OccupationNumber(NSites),1:mProj))

    n_avg_direct = mean(GFMC.get_obs_from_accumulator_bunching(ObsAccumulatorFile, 1))

    direct_mean = dropdims(mean( CObs.SaveConfigs,dims = (2,3)),dims=(2,3))
    
    @testset "ObservableAccumulator" begin

        @test isfile(outfile)
        GFMC.HDF5.h5open(outfile, "r") do file
            @test haskey(file, "OccupationNumber_denominator")

            OccupationNumber_denominator = read(file["OccupationNumber_denominator"])

            @test haskey(file, "OccupationNumber_numerator")
            OccupationNumber_numerator = read(file["OccupationNumber_numerator"])
            @test !iszero(OccupationNumber_denominator)
            @test !iszero(OccupationNumber_numerator)

            @test n_avg_direct[:,1] ≈ direct_mean atol = 1e-13
            @test n_avg ≈ n_avg_direct atol = 1e-10
        end
    end
end

@testitem "Accumulator TFI Discrete" begin
    include("utils.jl")
    using GreenFunctionMonteCarlo.Statistics:mean

    function E_critPoint_exact(L)
        return 1 - csc(pi / (2 * (2 * L + 1)))
    end
    σz(n::Bool) = (1 - 2 * n)
    σz(i, conf::AbstractArray) = σz(conf[i])

    using GreenFunctionMonteCarlo.LinearAlgebra
    NSites = 2
    NSteps = 1500
    mProj = 20
    NWalkers = 2
    nBranch = 5

    RNG = StableRNG(1234)

    Hilbert = BosonHilbertSpace(NSites, HardCoreConstraint())
    moves = eachcol(Bool.(I(NSites))) # each move flips a single spin
    offdiagElements = -ones(NSites)
    # Negative J: the diagonal (Ising) term has the opposite sign of the "Accumulator TFI" test above.
    # On this bipartite open chain, the staggered gauge transformation σz_i -> -σz_i on alternating sites
    # maps J -> -J while leaving the transverse-field term invariant, so the spectrum (and hence the exact
    # critical-point ground energy) is identical between positive- and negative-J TFI here.
    Hxx = DiagOperator(x-> -sum(σz(i, x) * σz(i + 1, x) for i in eachindex(x)[1:end-1]))

    H = localOperator(moves, offdiagElements, Hxx, Hilbert)

    config = BosonConfig(Hilbert)
    rand!(RNG,config)
    logψ = GFMC.EqualWeightSuperposition()

    # J is negative here, so H_xx(x) <= 0 for all reachable configurations, and Λ=1 keeps the self-loop
    # weight Λ - H_xx(x) safely positive without needing any pre-tuning of w_avg_estimate.
    DT = DiscretePropagator(1., nBranch)

    prob = GFMCProblem(config, NWalkers, DT; logψ, H, Hilbert)

    outfile = tempname()

    BObs = GFMC.BasicObserver(outfile, NSteps, NWalkers)
    CObs = GFMC.ConfigurationObserver(outfile, config, NSteps, NWalkers)

    BasicAccumulatorFile = GFMC.BasicAccumulator(outfile, mProj, NWalkers)
    ObsAccumulatorFile = GFMC.ObservableAccumulator(outfile,OccupationNumber(NSites),BasicAccumulatorFile, mProj, NWalkers, Threads.nthreads())

    Observer = GFMC.CombinedObserver((BObs, CObs,BasicAccumulatorFile, ObsAccumulatorFile))

    runGFMC!(prob, Observer, NSteps; rng = RNG)

    Energy = GFMC.getEnergies(BObs.TotalWeights, BObs.energies, mProj)
    Energy_direct = mean(GFMC.get_energy_from_accumulator_bunching(BasicAccumulatorFile, 1)) .* NSites
    @testset "BasicAccumulator" begin

        @test isfile(outfile)

        GFMC.HDF5.h5open(outfile, "r") do file
            @test haskey(file, "Gnp_denominator")

            Gnp_denominator = read(file["Gnp_denominator"])

            @test haskey(file, "en_numerator")
            en_numerator = read(file["en_numerator"])
            @test !iszero(Gnp_denominator)
            @test !iszero(en_numerator)

            @test Energy ≈ Energy_direct atol = 1e-10
            @test Energy[end÷2] ≈ E_critPoint_exact(NSites) rtol = 3e-2
        end
    end

    Gnps = GFMC.precomputeNormalizedAccWeight(BObs.TotalWeights,2mProj)

    n_avg = stack(getObs_diagonal(Gnps,CObs.SaveConfigs,BObs.reconfigurationTable,OccupationNumber(NSites),1:mProj))

    n_avg_direct = mean(GFMC.get_obs_from_accumulator_bunching(ObsAccumulatorFile, 1))

    direct_mean = dropdims(mean( CObs.SaveConfigs,dims = (2,3)),dims=(2,3))

    @testset "ObservableAccumulator" begin

        @test isfile(outfile)
        GFMC.HDF5.h5open(outfile, "r") do file
            @test haskey(file, "OccupationNumber_denominator")

            OccupationNumber_denominator = read(file["OccupationNumber_denominator"])

            @test haskey(file, "OccupationNumber_numerator")
            OccupationNumber_numerator = read(file["OccupationNumber_numerator"])
            @test !iszero(OccupationNumber_denominator)
            @test !iszero(OccupationNumber_numerator)

            @test n_avg_direct[:,1] ≈ direct_mean atol = 1e-13
            @test n_avg ≈ n_avg_direct atol = 1e-10
        end
    end
end

@testitem "RNG Threading reproducibility" begin
    include("utils.jl")
    using GreenFunctionMonteCarlo.Statistics:mean

    using GreenFunctionMonteCarlo.LinearAlgebra
    NSites = 8
    NSteps = 50
    mProj = 10
    NWalkers = 8

    Hilbert = BosonHilbertSpace(NSites, HardCoreConstraint())
    moves = eachcol(Bool.(I(NSites)))

    offdiagElements = -ones(NSites)
    Hxx = DiagOperator(x-> sum(σz(i, x) * σz(i + 1, x) for i in eachindex(x)[1:end-1]))

    H = localOperator(moves, offdiagElements, Hxx, Hilbert)

    config = BosonConfig(Hilbert)

    logψ = GFMC.EqualWeightSuperposition()
    CT = ContinuousTimePropagator(0.1)

    prob1 = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert)
    prob2 = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert)

    BObs1 = GFMC.BasicObserver(NSteps, NWalkers)
    BObs2 = GFMC.BasicObserver(NSteps, NWalkers)

    RNG1 = TaskLocalRNG()
    Random.seed!(RNG1, 1234)
    runGFMC!(prob1, BObs1, NSteps;rng = RNG1, parallelization = GFMC.MultiThreaded(2))

    RNG2 = TaskLocalRNG()
    Random.seed!(RNG2, 1234)
    runGFMC!(prob2, BObs2, NSteps;rng = RNG2, parallelization = GFMC.MultiThreaded(2))


    @testset "rng confs" begin
        @test prob1.Walkers.Configs == prob2.Walkers.Configs
        @test prob1.Walkers.WalkerWeights == prob2.Walkers.WalkerWeights
        @test BObs1.reconfigurationTable == BObs2.reconfigurationTable
        @test BObs1.energies == BObs2.energies
        @test BObs1.TotalWeights == BObs2.TotalWeights
    end
end

@testitem "MinimalReconfiguration records true walker ancestry" begin
    # Regression test for a bug where minimizeReconfiguration! computed the correct walker
    # ensemble (statistically) but left reconfigurationList (saved verbatim into
    # reconfigurationTable) out of sync with the copies actually performed by _replace_walkers!,
    # silently biasing every diagonal-observable backward-population trace (getPopulationMatrix!)
    # while leaving getEnergies (which never touches reconfigurationTable) unaffected.
    include("utils.jl")

    NWalkers = 20
    # getExampleHardcore uses ZeroDiagOperator() -- with EqualWeightSuperposition, every walker's
    # branching weight update is then identical regardless of its configuration, so
    # MinimalReconfiguration is (correctly) a documented no-op and this test would never exercise
    # any actual copies. getMinimalExample (a real TFI diagonal energy) gives genuine per-walker
    # weight spread instead.
    (;Hilbert,H) = getMinimalExample(6, 1.0, 1.0)

    config = BosonConfig(Hilbert)
    rand!(StableRNG(22), config)
    logψ = GFMC.EqualWeightSuperposition()
    CT = ContinuousTimePropagator(0.3)
    prob = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert, parallelization = GFMC.SingleThreaded())

    # Propagate a few real steps (no reconfiguration yet) so walker weights become genuinely
    # unequal -- MinimalReconfiguration is a documented no-op for identical weights, which
    # wouldn't exercise the copy/ancestry-recording logic this test targets at all.
    RNGs = [StableRNG(33)]
    for _ in 1:3
        GFMC.propagateWalkers!(prob.Walkers, H, logψ, Hilbert, CT, GFMC.SingleThreaded(), RNGs)
    end

    pre_configs = [copy(GFMC.getConfig(prob.Walkers, α)) for α in 1:NWalkers]

    reconfiguration = GFMC.MinimalReconfiguration(NWalkers)
    GFMC.reconfigurateWalkers!(prob.Walkers, reconfiguration, GFMC.SingleThreaded(), RNGs)
    # Sanity check that this scenario actually exercises real branching/death -- otherwise the
    # test below would trivially pass without ever touching the bug it's meant to catch.
    @test !isempty(reconfiguration.dead_walkers)

    list = GFMC.get_reconfigurationList(reconfiguration)
    @testset "reconfigurationList matches actual copies" begin
        for α in 1:NWalkers
            @test GFMC.getConfig(prob.Walkers, α) == pre_configs[list[α]]
        end
    end
end

@testitem "TFI L=3 open BC with Jastrow vs exact diagonalization" begin
    include("utils.jl")
    using GreenFunctionMonteCarlo.LinearAlgebra
    using GreenFunctionMonteCarlo.SparseArrays

    L = 3
    h, J = 0.4, 1.0

    Hilbert = BosonHilbertSpace(L, HardCoreConstraint())
    moves = [Bool[0 for _ in 1:L] for _ in 1:L]
    offdiagElements = zeros(Float64, L)
    for i in eachindex(moves)
        moves[i][i] = true
        offdiagElements[i] = -h
    end
    Hxx(conf) = -J * sum(σz(i, conf) * σz(i + 1, conf) for i in eachindex(conf)[1:end-1])
    H = localOperator(moves, offdiagElements, DiagOperator(Hxx), Hilbert)

    # Exact diagonalization reference (open BC; L=3 -> 2^3=8 states, trivial dense construction).
    # Convention: site 1 is the leftmost/most-significant tensor factor -- this only needs to be
    # self-consistent between the Hamiltonian and the Sz_i operators below, not matched to any
    # external convention.
    id = sparse(1.0I, 2, 2)
    σxmat = sparse([0.0 1.0; 1.0 0.0])
    σzmat = sparse([1.0 0.0; 0.0 -1.0])
    siteop(op, i) = foldl(kron, [k == i ? op : id for k in 1:L])
    Sz = [siteop(σzmat, i) for i in 1:L]
    Hmat = -J * sum(Sz[i] * Sz[i+1] for i in 1:(L-1)) - h * sum(siteop(σxmat, i) for i in 1:L)
    evals, evecs = eigen(Symmetric(Matrix(Hmat)))
    E0_exact = evals[1]
    ψ0 = evecs[:, 1]
    SzSz_exact = [ψ0' * (Sz[i] * Sz[j]) * ψ0 for i in 1:L, j in 1:L]

    # Decent (non-trivial) nearest-neighbor Jastrow guiding function on the hardcore-boson
    # representation -- EqualWeightSuperposition converges too slowly for a tight test.
    vij_jastrow = zeros(Float32, L, L)
    for i in 1:(L-1)
        vij_jastrow[i, i+1] = vij_jastrow[i+1, i] = 0.5f0
    end
    logψ = Jastrow(zeros(Float32, L), vij_jastrow)

    startConfig = BosonConfig(Hilbert)
    NWalkers = 300
    NSteps = 3000
    mProj = 30
    propagator = ContinuousTimePropagator(0.1)

    problem = GFMCProblem(startConfig, NWalkers, propagator; logψ, H, Hilbert)
    runGFMC!(problem, NoObserver(), 400; rng = StableRNG(101))
    ConfObs = ConfigObserver(startConfig, NSteps, NWalkers)
    runGFMC!(problem, ConfObs, NSteps; rng = StableRNG(202))

    E_gfmc = getEnergies(ConfObs, mProj)[end]
    @test isapprox(E_gfmc, E0_exact; atol = 0.1)

    SpinCorr = SpinCorrelations(L)
    SzSz_res = getObs_diagonal(ConfObs, SpinCorr, 1:mProj)
    SzSz_mat = get_matrix_from_tri(SzSz_res[end], SpinCorr.i_inds, SpinCorr.j_inds)
    @test all(isapprox.(SzSz_mat, SzSz_exact; atol = 0.05))
end

include("Jastrow_tests.jl")
include("parTemp_test.jl")
@run_package_tests verbose=true

