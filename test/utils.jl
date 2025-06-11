
using Random, StableRNGs, TestItems, GreenFunctionMonteCarlo, Test
using GreenFunctionMonteCarlo
import GreenFunctionMonteCarlo as GFMC
import GreenFunctionMonteCarlo.SmallCollections as SC

σz(n::Bool) = (1 - 2 * n)
σz(i, conf::AbstractArray) = σz(conf[i])

struct Hxx_TFI <: GFMC.DiagonalOperator
    J::Float64
end
(H::Hxx_TFI)(x::AbstractVector) = H.J* (sum(σz(i,x) * σz(i+1,x) for i in eachindex(x)[1:end-1]) - σz(x[begin]) * σz(x[end]))

function getMinimalExample(Nsites,h,J)
    Hilbert = BosonHilbertSpace(Nsites, HardCoreConstraint())

    moves = eachcol(GFMC.LinearAlgebra.I(Nsites))
    offdiag = -h .* ones(length(moves))

    diag_interaction = Hxx_TFI(J)

    H = localOperator(moves, offdiag, diag_interaction, Hilbert)
    
    return (; Hilbert, H)
end

get_move_type(::HardCoreConstraint) = Bool
get_move_type(C::OccupationNumberConstraint) = Int8

function getExample(Nsites,NMoves,rng,num_nonzero=nothing,constraint = HardCoreConstraint())
    Hilbert = BosonHilbertSpace(Nsites, constraint)

    getMove(num_nonzero::Nothing) = rand(rng,get_move_type(constraint),Nsites)
    getMove(num_nonzero::Int) = convert.(get_move_type(constraint,),SameLengthMove(Nsites,num_nonzero,rng))
    getMoves(NMoves,::HardCoreConstraint) = [getMove(num_nonzero) for i in 1:NMoves]
    function getMoves(NMoves,::OccupationNumberConstraint)
        moves = [getMove(num_nonzero) for i in 1:NMoves]
        append!(moves, -moves)
        return moves
    end
    moves = getMoves(NMoves,constraint)
    weights = -abs.(rand(rng,NMoves))
    if constraint isa OccupationNumberConstraint
        append!(weights, copy(weights))
    end
    H = localOperator(getMoves(NMoves,constraint), weights, ZeroDiagOperator(), Hilbert)
    
    return (; Hilbert, H)
end
getExampleHardcore(Nsites,NMoves,rng,num_nonzero=nothing) = getExample(Nsites,NMoves,rng,num_nonzero,HardCoreConstraint())

function SameLengthMove(Nsites, num_nonzero,rng=Random.default_rng())
    idxs = collect(1:Nsites)
    Random.shuffle!(rng,idxs)

    idx = idxs[1:num_nonzero]
    move = zeros(Bool,Nsites)
    move[idx] .= true
    return move
end 

function testSaveConf(SaveConfigs,TotalWeights,energies,reconfigurationTable,NSites,NWalkers,NSteps)
    @testset "SaveConfigs" begin

        @test size(SaveConfigs) == (NSites,NWalkers,NSteps)
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

function testPostMove(logψ,config,H,Hilbert)
    Buff = GFMC.allocate_GWF_buffer(logψ,config)

    move = H.moves[begin]
    for moveidx in eachindex(H.moves)
        move = H.moves[moveidx]
        GFMC.isapplicable(config,move,Hilbert) && break
    end
    xpr = copy(config)
    apply!(xpr,move)
    GFMC.post_move_affect!(Buff,xpr,move,logψ)

    Buff2 = GFMC.allocate_GWF_buffer(logψ,xpr)
    compareBuffers(Buff, Buff2)
end

compareBuffers(Buff::GFMC.SimpleJastrow_GWF_Buffer, Buff2::GFMC.SimpleJastrow_GWF_Buffer) = @test Buff.h_i ≈ Buff2.h_i  atol = 1e-14

function testRun(config,logψ,H,Hilbert;NWalkers = 20,NSteps = 10)

    CT = ContinuousTimePropagator(1.)

    RNG = StableRNG(1234)
    RNG_jastrow = StableRNG(1234)
    
    prob = GFMCProblem(config, NWalkers, CT; logψ=NaiveFunction(logψ), H, Hilbert,parallelization = GFMC.SingleThreaded())
    prob_jastrow = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert,parallelization = GFMC.SingleThreaded())
    
    ConfSaver = ConfigObserver(config, NSteps,NWalkers)
    ConfSaver_jastrow = ConfigObserver(config, NSteps,NWalkers)

    runGFMC!(prob, ConfSaver, NSteps;rng = RNG,logger = nothing)
    runGFMC!(prob_jastrow, ConfSaver_jastrow, NSteps;rng = RNG_jastrow,logger = nothing)

    (;SaveConfigs,TotalWeights,energies,reconfigurationTable) = GFMC.getObs(ConfSaver_jastrow)
    testSaveConf(SaveConfigs,TotalWeights,energies,reconfigurationTable,NSites,NWalkers,NSteps)

    Obs_reference = GFMC.getObs(ConfSaver)

    @test Obs_reference.SaveConfigs == SaveConfigs
    @test Obs_reference.TotalWeights ≈ TotalWeights atol = 1e-10
    @test Obs_reference.energies ≈ energies atol = 1e-10
    @test Obs_reference.reconfigurationTable == reconfigurationTable
end
