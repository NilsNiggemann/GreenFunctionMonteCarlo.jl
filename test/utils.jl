
using Random, StableRNGs, TestItems, GreenFunctionMonteCarlo, Test
import GreenFunctionMonteCarlo as GFMC
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