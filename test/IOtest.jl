@testitem "IOtest" begin
    using GreenFunctionMonteCarlo
    using HDF5
    params = Dict(
        "Omega" => 1,
        "delta" => 2,
        "misc" => Dict(:a => 1, :b => 2.0, :c => "test"), 
        "mean_TotalWeights" => 55.,
        "configs" => zeros(Int8,20,10),
    )
    outfile = tempname()
    @testset "save_param_dict" begin
        save_params_dict(outfile, params, mode="w")
    end
    @testset "load_param_dict" begin
        h5open(outfile, "r") do file
            read(file) == params
        end
    end
end