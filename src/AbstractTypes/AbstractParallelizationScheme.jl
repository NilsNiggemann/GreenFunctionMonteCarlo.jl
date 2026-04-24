abstract type AbstractParallelizationScheme end

struct SingleThreaded <: AbstractParallelizationScheme end
struct MultiThreaded <: AbstractParallelizationScheme
    nTasks::Int
end

num_tasks_default() = Threads.nthreads()
function num_tasks_default(NWalkers::Int)
    nthreads = Threads.nthreads()
    # return nthreads
    numTasks = max(1, round(Int, NWalkers / 10nthreads)) * nthreads
    return min(numTasks, NWalkers) # Ensure numTasks does not exceed NWalkers
end

"""
    BatchMultiThreaded{T} <: AbstractParallelizationScheme

Parallelization scheme that splits work into consecutive batches for execution
across multiple tasks. Uses Polyester for static task scheduling.

# Fields
- `nTasks::Int`: Number of tasks to schedule (should be a small multiple of the number of threads for optimal performance).
- `batches::Vector{T}`: Collection of iteration batches assigned to tasks. This should contain the full iteration range split into `nTasks` consecutive chunks.
"""
struct BatchMultiThreaded{T} <: AbstractParallelizationScheme
    nTasks::Int
    batches::Vector{T}
end
"""
    BatchMultiThreaded(nTasks::Int, NWalkers::Int)
Constructor for `BatchMultiThreaded` that takes the number of tasks and the corresponding batches.
- `nTasks::Int`: The number of tasks to schedule.
- `NWalkers::Int`: The total number of walkers to be processed, used to generate batches if not provided.
"""
function BatchMultiThreaded(nTasks::Int,NWalkers::Int)
    nTasks = max(1, nTasks)
    batches = collect(ChunkSplitters.chunks(1:NWalkers, n = nTasks, split= ChunkSplitters.Consecutive()))
    return BatchMultiThreaded(nTasks, batches)
end
