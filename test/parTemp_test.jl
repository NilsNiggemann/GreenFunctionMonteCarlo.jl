

@testitem "swap confs" begin
    include("utils.jl")
    
    NSites = 10

    (Hilbert1,H1) = getMinimalExample(NSites,1,-0.1) # slightly ferro
    logψ1 = Jastrow(NSites,Float64)

    (Hilbert2,H2) = getMinimalExample(NSites,1,0.1) # slightly anti ferro
    logψ2 = Jastrow(NSites,Float64)
    el¹ = GFMC.LocalEnergy(H1, logψ1, Hilbert1)
    el² = GFMC.LocalEnergy(H2, logψ2, Hilbert2)

    x1 = BosonConfig(Hilbert1)
    x1 .= [iseven(i) for i in 1:NSites]
    x2 = BosonConfig(Hilbert2) # x1 is more suitable for H2, since paramagnetic
    x2 .= 0

    P_ex = GFMC.log_replica_exchange_weight(el¹, el², x1, x2)
    P_ex_back = GFMC.log_replica_exchange_weight(el¹, el², x2, x1)

    @testset "swapping" begin
        @test Hilbert1 == Hilbert2

        @test el²(x1) < el¹(x1) # x1 is more suitable for H2 than for H1
        @test el¹(x2) < el²(x2)  # x2 is more suitable for H2 than for H1


        @test P_ex > 0 # exp(P_ex) > 1

        @test P_ex_back < 0 # exp(P_ex_back) < 1

        new_x1 = copy(x1)
        new_x2 = copy(x2)

        GFMC.inplaceSwap!(new_x1,new_x2)
        @test x1 == new_x2
        @test x2 == new_x1
    end

    NWalkers = 1
    rng = StableRNG(1234)

    dtau = 0.1
    P = ProblemEnsemble([GFMCProblem(x1, NWalkers, ContinuousTimePropagator(dtau); logψ = logψ1, H = H1, Hilbert = Hilbert1),
                GFMCProblem(x2, NWalkers, ContinuousTimePropagator(dtau); logψ = logψ2, H = H2, Hilbert = Hilbert2)])

    oldconfs = [copy(GFMC.getConfigs(p)[1]) for p in P.problems]


    swap_weights = GFMC.get_swap_weights(P)

    GFMC.replicaExchange!(P;rng)
    confs = [GFMC.getConfigs(p)[1] for p in P.problems]

    @testset "test_exchanges" begin
        @test swap_weights[1,2,1,1] == swap_weights[1,1,1,2]== exp(P_ex)
        @test swap_weights[1,1,1,1] == 1 # no swap

        @test oldconfs[1] == confs[2]
        @test oldconfs[1] == confs[2]
        @test oldconfs[2] == confs[1]
    end


end