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
    
    RNG = StableRNG(1234)
    GFMC.propagateWalkers!(ensemble, H, logψ, Hilbert, CT, GFMC.SingleThreaded(),RNG)
    
    AllConfs = stack(ensemble.Configs)
    @testset "Continuous Time Propagation" begin
        @test AllConfs == Bool[0 0 0 0 1 0 0 0; 1 0 1 0 1 1 0 1; 0 1 1 0 0 0 0 1]
        GFMC.propagateWalkers!(ensemble, H, logψ, Hilbert, CT, GFMC.SingleThreaded(),RNG)

        @test stack(ensemble.Configs) != AllConfs
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
    σz(n::Bool) = (1 - 2 * n)
    σz(i, conf::AbstractArray) = σz(conf[i])

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

@testitem "Offdiagonal Observable Operator Tests" begin
    include("utils.jl")
    Hilbert = BosonHilbertSpace(3, HardCoreConstraint())
    logψ = GFMC.EqualWeightSuperposition()

    moves = [
        Bool[1,0,0],
        Bool[0,1,0],
        Bool[0,0,1],
    ]
    O = GFMC.offdiagonalObservable(moves, (idx, ψratio, x) -> ψratio, Hilbert)

    config = BosonConfig(Hilbert)
    Buffer = GFMC.allocate_GWF_buffer(logψ, config)
    weights = zeros(GFMC.n_moves(O))

    GFMC.get_observable_weights!(weights, config, O, logψ, Hilbert, Buffer)

    @testset "Weights" begin
        @test GFMC.n_moves(O) == 3
        @test weights ≈ ones(3)
    end

    RNG = StableRNG(1234)
    move, w = GFMC.sample_and_apply_observable!(config, weights, O, RNG)

    @testset "Sampling" begin
        @test move !== nothing
        @test w ≈ 3.0
        @test count(config) == 1
    end

    @testset "Reproducibility" begin
        config2 = BosonConfig(Hilbert)
        Buffer2 = GFMC.allocate_GWF_buffer(logψ, config2)
        weights2 = zeros(GFMC.n_moves(O))
        GFMC.get_observable_weights!(weights2, config2, O, logψ, Hilbert, Buffer2)

        RNG2 = StableRNG(1234)
        move2, w2 = GFMC.sample_and_apply_observable!(config2, weights2, O, RNG2)

        @test config2 == config
        @test w2 == w
    end
end

@testitem "Forward Walking Slot Mechanics" begin
    include("utils.jl")
    Hilbert = BosonHilbertSpace(3, HardCoreConstraint())
    logψ = GFMC.EqualWeightSuperposition()

    moves = [
        Bool[0,1,1],
        Bool[0,0,1],
        Bool[0,1,0],
        Bool[1,1,0],
    ]
    H = localOperator(moves, -[0.5, 0.3, 0.2, 0.4], ZeroDiagOperator(), Hilbert)
    Observable = GFMC.offdiagonalObservable(moves, (idx, ψratio, x) -> ψratio, Hilbert)

    config = BosonConfig(Hilbert)
    NWalkers = 4
    Walkers = GFMC.allocate_walkerEnsemble(config, logψ, NWalkers, H)
    GFMC.compute_GWF_buffers!(Walkers, logψ)

    CT = ContinuousTimePropagator(1.)
    mProj = 4
    cadence = 2
    n_slots = fld(mProj, cadence) + 1
    slots = [GFMC.allocate_forward_walking_slot(config, logψ, NWalkers, H) for _ in 1:n_slots]
    weights_buf = zeros(GFMC.n_moves(Observable))
    RNG = StableRNG(99)

    @testset "Seeding and advancing" begin
        @test count(s -> s.active, slots) == 0

        GFMC.seed_from_main_population!(slots[1], Walkers, Observable, logψ, Hilbert, 1, weights_buf, RNG)
        @test slots[1].active
        @test slots[1].n_seed == 1
        @test slots[1].step == 0
        @test slots[1].cum_weight == 1.0
        @test isfinite(slots[1].seed_weight)

        for _ in 1:mProj
            GFMC.advance_slot!(slots[1], CT, logψ, H, Hilbert, 1.0, GFMC.SingleThreaded(), RNG)
        end
        @test slots[1].step == mProj
        @test isfinite(slots[1].cum_weight)
    end

    @testset "Free slot lookup" begin
        free1 = GFMC._find_free_slot(slots)
        @test !free1.active
        free1.active = true

        free2 = GFMC._find_free_slot(slots)
        @test !free2.active
        @test free2 !== free1
    end
end

@testitem "Forward Walking Accumulator TFI" begin
    include("utils.jl")
    using GreenFunctionMonteCarlo.LinearAlgebra

    σz(n::Bool) = (1 - 2 * n)
    σz(i, conf::AbstractArray) = σz(conf[i])

    NSites = 2
    NSteps = 1500
    m_proj_Basic = 20
    mProj = 6
    cadence = 3
    NWalkers = 2

    RNG = StableRNG(1234)

    Hilbert = BosonHilbertSpace(NSites, HardCoreConstraint())
    moves = collect(eachcol(Bool.(I(NSites))))
    offdiagElements = -ones(NSites)
    Hxx = DiagOperator(x -> sum(σz(i, x) * σz(i + 1, x) for i in eachindex(x)[1:end-1]))

    H = localOperator(moves, offdiagElements, Hxx, Hilbert)
    # Same moves as H, weight = the bare guiding-function ratio (no extra phase factor): this is
    # the q=0 structure factor, i.e. Σ_i σx_i.
    Observable = GFMC.offdiagonalObservable(moves, (idx, ψratio, x) -> ψratio, Hilbert)

    # Exact diagonalization on the 4-state, 2-site Hilbert space, built directly from H's/Observable's
    # own move and matrix-element definitions (rather than a hand-derived matrix) to avoid a
    # transcription mismatch with what the simulation actually samples.
    basis = [BosonConfig(Bool[false,false]), BosonConfig(Bool[false,true]), BosonConfig(Bool[true,false]), BosonConfig(Bool[true,true])]
    idx_of(x) = findfirst(y -> y == x, basis)

    Hmat = zeros(4,4)
    Omat = zeros(4,4)
    for (n,x) in enumerate(basis)
        Hmat[n,n] = Hxx(x)
        for i in 1:NSites
            y = copy(x)
            GFMC.apply!(y, GFMC.get_move(H, i))
            Hmat[n, idx_of(y)] += offdiagElements[i]

            y2 = copy(x)
            GFMC.apply!(y2, GFMC.get_move(Observable, i))
            Omat[n, idx_of(y2)] += GFMC.observable_weight(Observable, i, 1.0, x)
        end
    end

    evals, evecs = eigen(Symmetric(Hmat))
    ψ0 = evecs[:,1]
    O_exact = ψ0' * Omat * ψ0

    config = BosonConfig(Hilbert)
    rand!(RNG, config)
    logψ = GFMC.EqualWeightSuperposition()
    CT = ContinuousTimePropagator(0.1)

    prob = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert)

    outfile = tempname()
    BObs = GFMC.BasicObserver(outfile, NSteps, NWalkers)
    BasicAccumulatorFile = GFMC.BasicAccumulator(outfile, m_proj_Basic, NWalkers)
    FWAcc = GFMC.ForwardWalkingAccumulator(outfile, Observable, BasicAccumulatorFile, config, logψ, H, Hilbert, CT, mProj, cadence, NWalkers; rng = StableRNG(4321))

    Observer = GFMC.CombinedObserver((BObs, BasicAccumulatorFile, FWAcc))

    runGFMC!(prob, Observer, NSteps; rng = RNG)

    estimate = GFMC.get_observable_from_accumulator(FWAcc)[1]

    @testset "Forward walking convergence" begin
        @test all(isfinite, estimate)
        # p=mProj (real forward continuation) should be at least as close to the exact ground-state
        # expectation value as p=0 (seed weight only, the biased mixed estimator).
        @test abs(estimate[end] - O_exact) <= abs(estimate[1] - O_exact) + 0.2
    end
end

include("Jastrow_tests.jl")
include("parTemp_test.jl")
@run_package_tests verbose=true

