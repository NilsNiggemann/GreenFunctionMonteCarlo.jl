using GreenFunctionMonteCarlo
using Test
using Random
##
@testset "Bosonic Configuration Tests" begin
    Hilbert = BosonHilbertSpace(10, (OccupationNumberConstraint(0, 1),))

    config = BosonConfig(zeros(UInt8,10))
    
    @testset "UInt8 Variables" begin
        @test size(config) === (10,)
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

end