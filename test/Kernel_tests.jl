

@testitem "Kernels" begin
    include("utils.jl")
    (;Hilbert,H) = getExampleHardcore(3,4,StableRNG(1234))

    function convert_H(H,Nsites)
        numMoves = length(H.moves)
        moveMatrix = zeros(Bool, Nsites,numMoves)
        conf = BosonConfig(Hilbert)
        for (i, move) in enumerate(H.moves)
            newConf = copy(conf)
            GFMC.apply!(newConf, move)
            moveMatrix[:,i] .= newConf
        end
        return GFMC.MatrixMoveOperator(moveMatrix, H.off_diag, ZeroDiagOperator())
    end
    config = BosonConfig(Hilbert)
    Nsites = length(config)

    H = convert_H(H, Nsites)

    NWalkers = 10
    RNG = StableRNG(1234)
    rand!(RNG,config)

    vij = 0.1*rand(RNG,Nsites,Nsites)
    vij = vij + vij'

    jastrow(x,params) = GFMC.evaluate_jastrow(x,params[1],params[2])

    logψ = GFMC.ParametrizedFunction(jastrow,(zeros(Nsites), vij))

    WE = GFMC.allocate_walkerEnsemble(config, logψ, NWalkers, size(H.moves,2), GFMC.KernelParallelization(Array))

    @testset "ArrayWalkerEnsemble" begin

        @test WE isa GFMC.ArrayWalkerEnsemble
        @test isconcretetype(typeof(WE))
        @test all(x == config for x in eachcol(WE.Configs))
    end

    GFMC.get_markov_weights!(WE,H,logψ,Hilbert)
    @testset "markov Weights" begin

        @test !all(iszero,WE.MoveWeights)
        @test all(WE.MoveWeights .>= 0)
    end
end