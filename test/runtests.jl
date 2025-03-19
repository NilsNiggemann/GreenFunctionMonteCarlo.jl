using GreenFunctionMonteCarlo
import GreenFunctionMonteCarlo as GFMC
using Test
using Random
using StableRNGs
##

@testset "Bosonic Configuration Tests" begin
    Hilbert = BosonHilbertSpace(10, OccupationNumberConstraint(0, 1))

    config = BosonConfig(Hilbert)
    
    @testset "UInt8 Variables" begin
        @test size(config) === (10,)
        @test config == BosonConfig(zeros(UInt8,10))
        @test fulfills_constraint(config, Hilbert)
        rand!(config)
        @test !fulfills_constraint(config, Hilbert)
    end

    config1 = BosonConfig(zeros(Bool,10))
    config2 = BosonConfig(BitVector(zeros(Bool,10)))
    
    for config in (config1, config2)
        @testset "Bool Variables" begin
            @test size(config) === (10,)
            @test fulfills_constraint(config, Hilbert)
            rand!(config)
            @test fulfills_constraint(config, Hilbert)
        end
    end

    @testset "HardCoreConstraint" begin
        Hilbert = BosonHilbertSpace(10, HardCoreConstraint())
        config = BosonConfig(Hilbert)
        @test fulfills_constraint(config, Hilbert)
        rand!(config)
        @test fulfills_constraint(config, Hilbert)
    end
end
##


@testset "Sparse Move Tests" begin
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

    config = BosonConfig(UInt8[1,1,0])
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

@testset "Bosonic Walker Ensemble Tests" begin
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
@testset "main Usage" begin
    
    Hilbert = BosonHilbertSpace(3, HardCoreConstraint())
    config = BosonConfig(Hilbert)
    rand!(config)
    logψ = GFMC.EqualWeightSuperposition()
    NWalkers = 10
    CT = ContinuousTimePropagator(1.)
    
    RNG = StableRNG(1234)
    H = localOperator(
        [
            Bool[0,1,1],
            Bool[0,0,1],
            Bool[0,1,0],
            Bool[1,1,0],
            ]
    , -[0.5, 0.3, 0.2 ,0.4], ZeroDiagOperator(), Hilbert)

    prob = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert)

    @testset "runGFMC" begin
        obs = NoObserver()
        runGFMC!(prob, obs, 10, RNG)
        AllConfs = stack(prob.WE.Configs)
        @test AllConfs != zeros(Bool,3,NWalkers)
    end

    outfile = tempname()

    NSteps = 5

    ConfSaverFile = GFMC.configObserver(outfile, config, NSteps,NWalkers)

    @testset "ConfigObserver" begin

        @test isfile(outfile)
        runGFMC!(prob, ConfSaverFile, NSteps, RNG)

        GFMC.HDF5.h5open(outfile, "r") do file
            println(file)
            println(keys(file))
            @testset "SaveConfigs" begin
                @test haskey(file, "SaveConfigs")

                SaveConfigs = read(file["SaveConfigs"])
                @test size(SaveConfigs) == (3,NWalkers,NSteps)
                @test eltype(SaveConfigs) == Bool
                @test !iszero(SaveConfigs)
            end
            @testset "TotalWeights" begin
                @test haskey(file, "TotalWeights")
                TotalWeights = read(file["TotalWeights"])
                @test size(TotalWeights) == (NSteps,)
                @test eltype(TotalWeights) == Float64
                @test !iszero(TotalWeights)
            end
            @testset "energies" begin
                @test haskey(file, "energies")
                energies = read(file["energies"])
                @test size(energies) == (NSteps,)
                @test eltype(energies) == Float64
                @test !iszero(energies)
            end
            @testset "reconfigurationTable" begin
                @test haskey(file, "reconfigurationTable")
                reconfigurationTable = read(file["reconfigurationTable"])
                @test size(reconfigurationTable) == (NWalkers,NSteps)
                @test eltype(reconfigurationTable) == Int
                @test !iszero(reconfigurationTable)
            end
        end
    end
end