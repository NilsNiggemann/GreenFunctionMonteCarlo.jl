abstract type AbstractParallelizationScheme end

struct SingleThreaded <: AbstractParallelizationScheme end
struct MultiThreaded <: AbstractParallelizationScheme
    nTasks::Int
end

num_tasks_default() = Threads.nthreads()
function num_tasks_default(NWalkers::Int)
    nthreads = Threads.nthreads()
    return nthreads
    numTasks = max(1, round(Int, NWalkers / 10nthreads)) * nthreads
    return min(numTasks, NWalkers) # Ensure numTasks does not exceed NWalkers
end
