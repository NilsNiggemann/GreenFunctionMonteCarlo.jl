abstract type AbstractParallelizationScheme end

struct SingleThreaded <: AbstractParallelizationScheme end
struct MultiThreaded <: AbstractParallelizationScheme
    nTasks::Int
end

num_tasks_default() = Threads.nthreads()
function num_tasks_default(NWalkers::Int)
    nthreads = Threads.nthreads()
    numTasks = max(1, round(Int, 2* NWalkers / nthreads)) * nthreads
    return min(numTasks, NWalkers) # Ensure numTasks does not exceed NWalkers
end
