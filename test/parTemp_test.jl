

@testitem "swap confs" begin
    include("utils.jl")
    
    NSites = 10

    function coup(J,i,j)
        i,j = sort(GFMC.SA.SVector(i,j))
        i == j-1 && return 0.5J
        i == 1 && j == Nsites && return 0.5J
        return 0.
    end

    (Hilbert1,H1) = getMinimalExample(NSites,1,-0.1) # slightly ferro
    logψ1 = Jastrow(NSites,Float64)
    # logψ1.v_ij .= -0.1*[coup(H1.diag.J,i,j) for i in 1:NSites, j in 1:NSites]

    (Hilbert2,H2) = getMinimalExample(NSites,1,0.1) # slightly anti ferro
    logψ2 = Jastrow(NSites,Float64)
    # logψ2.v_ij .= -0.1*[coup(H1.diag.J,i,j) for i in 1:NSites, j in 1:NSites]
    el¹ = GFMC.LocalEnergy(H1, logψ1, Hilbert1)
    el² = GFMC.LocalEnergy(H2, logψ2, Hilbert2)

    x1 = BosonConfig(Hilbert1)
    x1 .= [iseven(i) for i in 1:NSites]
    x2 = BosonConfig(Hilbert2) # x1 is more suitable for H2, since paramagnetic

    x2 .= 0

    @testset "swapping" begin
        @test Hilbert1 == Hilbert2

        @test el²(x1) < el¹(x1) # x1 is more suitable for H2 than for H1
        @test el¹(x2) < el²(x2)  # x2 is more suitable for H2 than for H1
        P_ex = GFMC.P_exchange(el¹, el², x1, x2)

        @test P_ex == 1.0

        P_ex_back = GFMC.P_exchange(el¹, el², x2, x1)
        @test P_ex_back < P_ex

        println("P_ex = ", P_ex)
    end

end