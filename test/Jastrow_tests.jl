@testitem "HardcoreConstraint, fixed move length" begin
    include("utils.jl")
    
    NSites = 10
    RNG = StableRNG(1234)
    logψ = Jastrow(NSites,Float64)
    params = get_params(logψ)
    rand!(RNG,params)
    logψ.v_ij .= GFMC.LinearAlgebra.Symmetric(logψ.v_ij)
    
    (Hilbert,H) = getExample(NSites,4,RNG,3,HardCoreConstraint())

    config = BosonConfig(Hilbert)
    rand!(RNG,config)

    @testset "Wavefunction Ratio" begin
        TestWFRatio(logψ,config,H,Hilbert)
    end
    
    @testset "post_move" begin
        testPostMove(logψ,config,H,Hilbert)
    end
    
    @testset "run Jastrow" begin
        testRun(config,logψ,H,Hilbert;NWalkers = 20,NSteps = 10)
    end
end

@testitem "HardcoreConstraint, variable move length" begin
    include("utils.jl")
    
    NSites = 10
    RNG = StableRNG(1234)
    logψ = Jastrow(NSites,Float64)
    params = get_params(logψ)
    rand!(RNG,params)
    logψ.v_ij .= GFMC.LinearAlgebra.Symmetric(logψ.v_ij)
    
    (Hilbert,H) = getExample(NSites,4,RNG,nothing,HardCoreConstraint())

    config = BosonConfig(Hilbert)
    rand!(RNG,config)

    @testset "Wavefunction Ratio" begin
        TestWFRatio(logψ,config,H,Hilbert)
    end
    
    @testset "post_move" begin
        testPostMove(logψ,config,H,Hilbert)
    end
    
    @testset "run Jastrow" begin
        testRun(config,logψ,H,Hilbert;NWalkers = 20,NSteps = 10)
    end
end


@testitem "OccupationNumberConstraint, fixed move length" begin
    include("utils.jl")
    
    NSites = 10
    RNG = StableRNG(1234)
    logψ = Jastrow(NSites,Float64)
    params = get_params(logψ)
    rand!(RNG,params)
    logψ.v_ij .= GFMC.LinearAlgebra.Symmetric(logψ.v_ij)
    
    (Hilbert,H) = getExample(NSites,4,RNG,3,OccupationNumberConstraint(0,1))

    config = BosonConfig(Hilbert)
    rand!(RNG,config,0:1)

    @testset "Wavefunction Ratio" begin
        TestWFRatio(logψ,config,H,Hilbert)
    end
    
    @testset "post_move" begin
        testPostMove(logψ,config,H,Hilbert)
    end
    
    @testset "run Jastrow" begin
        testRun(config,logψ,H,Hilbert;NWalkers = 20,NSteps = 10)
    end
end


@testitem "OccupationNumberConstraint, variable move length" begin
    include("utils.jl")
    
    NSites = 10
    RNG = StableRNG(1234)
    logψ = Jastrow(NSites,Float64)
    params = get_params(logψ)
    
    rand!(RNG,params)
    logψ.v_ij .= GFMC.LinearAlgebra.Symmetric(logψ.v_ij)
    
    (Hilbert,H) = getExample(NSites,4,RNG,3,OccupationNumberConstraint(0,1))

    config = BosonConfig(Hilbert)
    rand!(RNG,config,0:1)

    @testset "Wavefunction Ratio" begin
        TestWFRatio(logψ,config,H,Hilbert)
    end
    
    @testset "post_move" begin
        testPostMove(logψ,config,H,Hilbert)
    end
    
    @testset "run Jastrow" begin
        testRun(config,logψ,H,Hilbert;NWalkers = 20,NSteps = 10)
    end
end
