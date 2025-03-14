using GreenFunctionMonteCarlo
import GreenFunctionMonteCarlo as GFMC
using Test
using Random
using SparseArrays
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

##

@testset "Sparse Move Tests" begin
    moves = sparse(Int8[
        0 -1 1;
        0 1 -1;
        1 0 1;
        -1 0 -1;
    ]')
    moves_binary = sparse(Bool[
        0 1 1;
        1 0 1;
    ]')

    weights = [0.5, 0.5, 0.2,0.2]
    weights_bin = [0.5, 0.2]

    diag = ZeroDiagOperator()
    local_op = localOperator(moves, weights, diag)

    config = BosonConfig(UInt8[1,1,0])
    config_Bin = BosonConfig(Bool[1,1,0])

    @testset "SparseMove" begin
        move = GFMC.get_move(local_op, 1)
        GFMC.apply!(config, move)
        @test config == [1,0,1]
    end

    @testset "SparseMoveBinary" begin
        local_op_binary = localOperator(moves_binary, weights_bin, diag)
        move_binary = GFMC.get_move(local_op_binary, 1)
        GFMC.apply!(config_Bin, move_binary)
        @test config_Bin == Bool[1,0,1]
    end
end