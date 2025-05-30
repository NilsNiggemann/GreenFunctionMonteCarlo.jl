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
    
    rand!(x)

    RNG = StableRNG(1234)
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
    rand!(config)
    logψ = GFMC.EqualWeightSuperposition()
    NWalkers = 10
    CT = ContinuousTimePropagator(1.)
    
    RNG = StableRNG(1234)

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

@testitem "BasicAccumulator TFI" begin
    include("utils.jl")

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

    outfile = tempname()

    BasicAccumulatorFile = GFMC.BasicAccumulator(outfile, mProj, NWalkers)
    BObs = GFMC.BasicObserver(outfile, NSteps, NWalkers)

    Observer = GFMC.CombinedObserver((BObs, BasicAccumulatorFile))
    runGFMC!(prob, NoObserver(), 200; rng = RNG)
    runGFMC!(prob, Observer, NSteps; rng = RNG)

    Energy = GFMC.getEnergies(BObs.TotalWeights, BObs.energies, mProj)
    Energy_direct = BasicAccumulatorFile.en_numerator ./ BasicAccumulatorFile.Gnp_denominator .*NSites

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

            @test Energy[end÷2] ≈ E_critPoint_exact(NSites) rtol = 2e-2
        end
    end
end


include("Jastrow_tests.jl")
@run_package_tests verbose=true