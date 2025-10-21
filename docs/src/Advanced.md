# Advanced Usage
## Accumulators
For long- running simulations with many walkers, the amount of data generated can be substantial. The `ConfigSaver`, introduced in [Transverse Field Ising Model Example](Example_transverseFieldIsing.md), must store the configurations of *each* walker at *each* observation step, which can lead to high memory consumption. Typically, one is ultimately interested in the expectation values of a handful of observables, which may also be obtained by accumulating data on-the-fly during the simulation.
GreenFunctionMonteCarlo.jl provides a [`BasicAccumulator`](@ref) which stores the energy and some additional data required for other accumulators, as well as an  [`ObservableAccumulator`](@ref) which can be used to accumulate arbitrary observables defined as subtypes of `AbstractObservable`.

In summary, the advantages of accumulators are:
- **Reduced Storage Usage**
The drawbacks are:
- **Limited Post-Processing Flexibility**: The observables cannot be obtained from post-processing. The maximum projection time needs to be defined before the simulation.
- **Increased memory usage**: To avoid costly recomputation of expensive obersvables for the same configurations, accumulators use buffers to store the observables for the last couple of iterations to do the projection. The memory size required is given by `NObs* NWalkers * N_projection_steps*size(data_type)`. It can be advisable to use a more efficient `data_type` for the buffer, such as `Float32`. For example, accumulating the local magnetization for *N_sites*, the memory usage will be
```@example julia
N_sites = 500
N_obs = N_sites
NWalkers = 5000
m_proj = 150
Base.format_bytes(Base.summarysize(zeros(Float32, N_obs))* NWalkers * m_proj)
```
- **Possible Numerical Instabilities**: Accumulating data on-the-fly can introduce numerical instabilities, especially when the estimated average weight of the walkers is not accurate (in which case at each step, some very large and very small numbers are added to numerator and denominator). It is advisable to initialize the `ContinuousTimePropagator` with a reasonable estimate of the ground state energy to mitigate this issue. Below, it is shown how to do that. Another strategy is to use more than one bin for the accumulation, so that not too many numbers are added up to the same storage.

### Example run
We first run a simulation without measuring observables to estimate the average weight of the walkers. It can be used as a replacement for equilibrating the walkers with a first `runGFMC!` call. With each epoch, the estimate will be a bit better, but typically not very many should be required.
```julia
problem = GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ, H, Hilbert)
mean_TotalWeights, w_avg_estimate = GFMC.estimate_weights_continuousTime!(problem;Nepochs=4,Nsamples=400,verbose=true,logger = nothing)
```

The `BasicAccumulator` stores running sums of the numerator and the denominator for the energy. It also provides internal buffers which are used for the reconfiguration process.
In order to measure other observables, one can use the `ObservableAccumulator`, which holds a reference to the `BasicAccumulator` to access the internal buffers. Here, we measure the local magnetization and the density-density correlation function as an example.

The set of observables is packed as a tuple into a `CombinedObserver`, which makes sure that each observable is updated at each observation step.

```julia
outfile = "observables.h5"
num_bins = 32 # number of bins for the accumulation
BAcc = BasicAccumulator(outfile,m_proj, NWalkers;weight_normalization = mean_TotalWeights, num_bins = num_bins , bin_elements = NSteps_total ÷ num_bins)

OccNum_Acc = ObservableAccumulator(outfile,OccupationNumber(Nsites),BAcc, mProj, NWalkers, Threads.nthreads())

Observer = CombinedObserver((BAcc, OccNum_Acc, )) # note the tuple. Include more observable accumulators as needed

CT = ContinuousTimePropagator(dtau, w_avg_estimate) # use the estimated average weight here
runGFMC!(problem, Observer, NSteps_total; Propagator = CT) # we may override the propagator from "problem" by passing a different one as a keyword argument
```

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

    Obs = (;Energy, OccupationNumbers)
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
save_params_dict(outfile, params, mode="w")
h5open(outfile, "r") do f
    println(read(f))
end
```