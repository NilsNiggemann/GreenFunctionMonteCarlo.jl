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
"""
struct BatchMultiThreaded <: AbstractParallelizationScheme
    nTasks::Int
end