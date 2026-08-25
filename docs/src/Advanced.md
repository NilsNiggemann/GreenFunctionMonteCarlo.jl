# Advanced Usage
## Accumulators
For long- running simulations with many walkers, the amount of data generated can be substantial. The `ConfigSaver`, introduced in [Transverse Field Ising Model Example](Example_transverseFieldIsing.md), must store the configurations of *each* walker at *each* observation step, which can lead to high memory consumption. Typically, one is ultimately interested in the expectation values of a handful of observables, which may also be obtained by accumulating data on-the-fly during the simulation.
GreenFunctionMonteCarlo.jl provides a [`BasicAccumulator`](@ref) which stores the energy and some additional data required for other accumulators, as well as a [`LazyObservableAccumulator`](@ref) and an [`ObservableAccumulator`](@ref) which can be used to accumulate arbitrary observables defined as subtypes of `AbstractObservable`.

In summary, the advantages of accumulators are:
- **Reduced Storage Usage**
The drawbacks are:
- **Limited Post-Processing Flexibility**: The observables cannot be obtained from post-processing. The maximum projection time needs to be defined before the simulation.
- **Possible Numerical Instabilities**: Accumulating data on-the-fly can introduce numerical instabilities, especially when the estimated average weight of the walkers is not accurate (in which case at each step, some very large and very small numbers are added to numerator and denominator). It is advisable to initialize the `ContinuousTimePropagator` with a reasonable estimate of the ground state energy to mitigate this issue. Below, it is shown how to do that. Another strategy is to use more than one bin for the accumulation, so that not too many numbers are added up to the same storage.

### Example run
A full worked example — estimating the average walker weight, setting up a `BasicAccumulator` together with `LazyObservableAccumulator`s for [`OccupationNumber`](@ref) and [`SpinCorrelations`](@ref), and running the simulation — is shown in the ["Efficient observable accumulation for large runs"](Example_transverseFieldIsing.md#Efficient-observable-accumulation-for-large-runs) section of the tutorial. The `BasicAccumulator` stores running sums of the numerator and the denominator for the energy, and also provides internal buffers used by `LazyObservableAccumulator`s and by the reconfiguration process; several observables are combined via a `CombinedObserver`, which updates each of them at every observation step.

In the tutorial's example, `outfile = "observables.h5"`, `m_proj = mProj`, and `NSteps_total = 2000`.

### Example evaluation
After the simulation, the accumulated data can be accessed from the HDF5 file. We still need to do some post-processing, such as averaging over the bins. This is not done automatically in the previous step, because the binning may also be used to estimate error bars, provided that the simulation is long enough to decorrelate the bins.

```julia
# read the data from the HDF5 file and pack the numerators and denominators into tuples that mock the structure of the accumulators.
# Note: If the accumulators themselves are still in memory, we can of course also use them instead.
function readObs(file,bunching)
    BasicAccMock = (;en_numerator = h5read(file, "en_numerator"),Gnp_denominator = h5read(file, "Gnp_denominator"))

    Energy = GFMC.get_energy_from_accumulator_bunching(BasicAccMock, bunching)

    OccNumMock = (;Obs_numerators = h5read(file, "OccupationNumber_numerator"),
    Obs_denominators = h5read(file, "OccupationNumber_denominator"))
    OccupationNumbers = GFMC.get_obs_from_accumulator_bunching(OccNumMock, bunching)
    m_values = h5read(file, "OccupationNumber_m_values")
    Obs = (;Energy, OccupationNumbers, m_values)
    return Obs
end
```
If we call the function with `bunching = 8`, the following will happen:
- Numerators and denominators will be normalized to avoid numerical instabilities when adding up many numbers.
- Numerators and denominators of each 8 bins will be added up.
- We will divide the summed numerators by the summed denominators to obtain the final expectation values.
- Since there were `num_bins = 32` bins in total, we will obtain 4 data points for each observable, which can be used to estimate error bars.

!!! warning
    To estimate error bars, you must ensure that the bins are sufficiently decorrelated. This typically requires long simulations.

## Saving simulation parameters
To store parameters of the simulation, such as model parameters, number of walkers, time step, etc., GreenFunctionMonteCarlo.jl provides a function `save_params_dict`, which saves a dictionary of key-value pairs to an HDF5 file. This can be useful for keeping track of the simulation settings alongside the results.
Nested dicts are also supported and will be stored as HDF5 groups.
```@example julia
using GreenFunctionMonteCarlo
using GreenFunctionMonteCarlo.HDF5
params = Dict(
    "mProj" => 5,
    "misc" => Dict(:a => 1, :b => 2.0, :c => "test"), 
    "NSteps_total" => 55,
    "mean_TotalWeights" => 55.,
    "configs" => zeros(Int8,3,2),
)
outfile = tempname()
save_params_dict(outfile, Dict("params" => params), mode="w")
h5open(outfile, "r") do f
    println(read(f))
end
```

## Buffered ObservableAccumulator
If the observables are costly to compute, it can be advisable to use the [`ObservableAccumulator`](@ref) instead of [`LazyObservableAccumulator`](@ref).
- **Increased memory usage**: To avoid costly recomputation of expensive obersvables for the same configurations, accumulators use buffers to store the observables for the last couple of iterations to do the projection. The memory size required is given by `NObs* NWalkers * N_projection_steps*size(data_type)`. It can be advisable to use a more efficient `data_type` for the buffer, such as `Float32`. For example, accumulating the local magnetization for *N_sites*, the memory usage will be
```@example julia
N_sites = 500
N_obs = N_sites
NWalkers = 5000
m_proj = 150
Base.format_bytes(Base.summarysize(zeros(Float32, N_obs))* NWalkers * m_proj)
```

## Parallelization
`GFMCProblem` accepts a `parallelization` keyword argument controlling how walkers are distributed across threads, defaulting to `MultiThreaded` (which uses `Threads.@spawn` over chunks of walkers). Two other schemes are available:
- `SingleThreaded`: no parallelization, mainly useful for debugging.
- [`BatchMultiThreaded`](@ref): uses [Polyester.jl](https://github.com/JuliaSIMD/Polyester.jl) for static, low-overhead task scheduling instead of `Threads.@spawn`. It generally performs better than `MultiThreaded` when many CPU cores are available, and is used throughout the tutorial's accumulator example.

```julia
parallelization = BatchMultiThreaded(Threads.nthreads()) # nTasks should be a small multiple of Threads.nthreads()
problem = GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ, H, Hilbert, parallelization)
```

The parallelization scheme can also be overridden per call, e.g. `runGFMC!(problem, Observer, NSteps; parallelization = BatchMultiThreaded(Threads.nthreads()))`.

## DiscretePropagator
[`ContinuousTimePropagator`](@ref) is the recommended propagator for most problems. [`DiscretePropagator`](@ref) implements the more traditional discrete-time projector `G = Λ·I - H`: at each of `nBranch` sub-steps, a walker either stays put (self-loop weight `Λ - H_xx(x)`) or moves via an off-diagonal matrix element. This requires choosing `Λ ≥ H_xx(x)` for every configuration visited, and is only efficient when a reasonably tight such `Λ` is known — i.e. when the diagonal part `H_xx` is small or zero. Otherwise, the self-loop probability dominates the sampling and `ContinuousTimePropagator` should be preferred.

```julia
propagator = DiscretePropagator(Λ, nBranch) # Λ ≥ max(Hxx(x)) over all reachable configurations x
```
