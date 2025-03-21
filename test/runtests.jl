using GreenFunctionMonteCarlo
import GreenFunctionMonteCarlo as GFMC
using Test
using Random
using StableRNGs
import GreenFunctionMonteCarlo.SmallCollections as SC

function getExampleHardcore(Nsites,NMoves,rng,num_nonzero=nothing)
    Hilbert = BosonHilbertSpace(Nsites, HardCoreConstraint())

    getMove(num_nonzero::Nothing) = rand(rng,Bool,Nsites)
    getMove(num_nonzero::Int) = SameLengthMove(Nsites,num_nonzero,rng)
    H = localOperator(
        [
            getMove(num_nonzero) for i in 1:NMoves
        ]
    , -abs.(rand(NMoves)), ZeroDiagOperator(), Hilbert)
    
    return (; Hilbert, H)
end

function SameLengthMove(Nsites, num_nonzero,rng=Random.default_rng())
    idxs = collect(1:Nsites)
    Random.shuffle!(rng,idxs)

    idx = idxs[1:num_nonzero]
    move = zeros(Bool,Nsites)
    move[idx] .= true
    return move
end

##

@testset "Bosonic Configuration Tests" begin
    Hilbert = BosonHilbertSpace(10, OccupationNumberConstraint(0, 1))
    RNG = StableRNG(1234)
    config = BosonConfig(Hilbert)
    
    @testset "UInt8 Variables" begin
        @test size(config) === (10,)
        @test config == BosonConfig(zeros(UInt8,10))
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
    @testset "FlipMoves" begin
        move = H.moves[1]
        @test move === GFMC.FlipMove(SC.SmallVector{2,Int}((2,3)))

        GFMC.move_dx(H.moves[1],config) === SC.SmallVector{2,Bool}((1,1))
        @test iszero(config)

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
function testSaveConf(SaveConfigs,TotalWeights,energies,reconfigurationTable,NSites,NWalkers,NSteps)
    @testset "SaveConfigs" begin

        @test size(SaveConfigs) == (NSites,NWalkers,NSteps)
        @test eltype(SaveConfigs) == Bool
        @test !iszero(SaveConfigs)
    end
    @testset "TotalWeights" begin
        @test size(TotalWeights) == (NSteps,)
        @test eltype(TotalWeights) == Float64
        @test !iszero(TotalWeights)
    end
    @testset "energies" begin
        @test size(energies) == (NSteps,)
        @test eltype(energies) == Float64
        @test !iszero(energies)
    end
    @testset "reconfigurationTable" begin
        @test size(reconfigurationTable) == (NWalkers,NSteps)
        @test eltype(reconfigurationTable) == Int
        @test !iszero(reconfigurationTable)
    end
end


@testset "main Usage" begin
    
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
        runGFMC!(prob, obs, 10, RNG)
        AllConfs = stack(prob.WE.Configs)
        @test AllConfs != zeros(Bool,3,NWalkers)
    end

    outfile = tempname()

    NSteps = 5

    ConfSaverFile = ConfigObserver(outfile, config, NSteps, NWalkers)

    @testset "ConfigObserver" begin

        @test isfile(outfile)
        runGFMC!(prob, ConfSaverFile, NSteps, RNG)

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

##

function TestWFRatio(logψ,conf,H,Hilbert;tol=1e-10)
    Buff = GFMC.allocate_GWF_buffer(logψ,conf)
    psiname = GFMC.guidingfunc_name(logψ)

    GFMC.compute_GWF_buffer!(Buff,logψ,conf)

    @testset "$psiname ψ(x´) / ψ(x)" begin
        for m in H.moves
            GFMC.isapplicable(conf,m,Hilbert) || continue
            psidiff = GFMC.log_psi_diff(conf, m, logψ, Buff, Hilbert)
            
            xpr = copy(conf)
            apply!(xpr,m)
            
            @test psidiff ≈ logψ(xpr) - logψ(conf) atol=tol
            
        end
    end
end

##
@testset "Jastrow Tests" begin
    
    @testset for num_nonzero in (nothing,3)
        RNG = StableRNG(1234)
        NSites = 20
        (Hilbert,H) = getExampleHardcore(NSites,4,RNG,num_nonzero)

        config = BosonConfig(Hilbert)
        
        rand!(RNG,config)

        logψ = Jastrow(config,Float64)

        params = get_params(logψ)
        rand!(RNG,params)
        logψ.v_ij .= GFMC.LinearAlgebra.Symmetric(logψ.v_ij)

        @testset "Wavefunction Ratio" begin
            TestWFRatio(logψ,config,H,Hilbert)
        end
        
        @testset "run Jastrow" begin
            NWalkers = 20
            NSteps = 10
            CT = ContinuousTimePropagator(1.)
            prob = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert,parallelization = GFMC.SingleThreaded())
        
            ConfSaver = ConfigObserver(config, NSteps,NWalkers)
        
            runGFMC!(prob, ConfSaver, NSteps, RNG)

            (;SaveConfigs,TotalWeights,energies,reconfigurationTable) = ConfSaver
            testSaveConf(SaveConfigs,TotalWeights,energies,reconfigurationTable,NSites,NWalkers,NSteps)
        end
    end
end