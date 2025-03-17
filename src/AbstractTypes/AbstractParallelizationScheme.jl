abstract type AbstractParallelizationScheme end

struct MultiThreaded <: AbstractParallelizationScheme
    nWorkChunks::Int
end