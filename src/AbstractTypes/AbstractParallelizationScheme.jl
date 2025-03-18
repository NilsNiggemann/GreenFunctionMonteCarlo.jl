abstract type AbstractParallelizationScheme end

struct SingleThreaded <: AbstractParallelizationScheme end
struct MultiThreaded <: AbstractParallelizationScheme
    nWorkChunks::Int
end
