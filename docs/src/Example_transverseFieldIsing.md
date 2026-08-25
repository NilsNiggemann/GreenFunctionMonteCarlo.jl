# Example: 1D Transverse Field Ising Model

## Defining the Hamiltonian
The transverse field Ising model is a quantum spin model that exhibits a quantum phase transition. The Hamiltonian for the model is given by:
```math
H = -J \sum_i^L \sigma_i^z \sigma_{i+1}^z - h \sum^L_i \sigma_i^x
```

where:
-  $L$ is the number of spins.
-  $J$ is the interaction strength between neighboring spins.
-  $h$ is the transverse field strength.
-  $\sigma^z$ and $\sigma^x$ are Pauli matrices.
At the critical point $h=1$, the energy for open boundary conditions is given by:
```math
E_{crit} = 1 - \cosec \left(\frac{\pi}{2(2L+1)}\right)
```
Below is an example of how to simulate the transverse field Ising model using `GreenFunctionMonteCarlo.jl`:

First, we implement our Hamiltonian. Since our problem is equivalent to a hardcore boson problem, we represent our spin configurations by Booleans for optimal performance.
```@example TFI
using GreenFunctionMonteCarlo

σz(n::Bool) = (1 - 2 * n)
σz(i, conf::AbstractArray) = σz(conf[i])

function transverse_field_ising(Nsites, h, J; periodic = false)
    Hilbert = BosonHilbertSpace(Nsites, HardCoreConstraint())

    moves = [Bool[0 for _ in 1:Nsites] for _ in 1:Nsites]
    offdiagElements = zeros(Float64, Nsites)

    for i in eachindex(moves)
        moves[i][i] = true
        offdiagElements[i] = -h
    end

    function Hxx(conf)
        E = -J * sum(σz(i, conf) * σz(i + 1, conf) for i in eachindex(conf)[1:end-1])
        if periodic
            E += -J * σz(Nsites, conf) * σz(1, conf)
        end
        return E
    end

    H = localOperator(moves, offdiagElements, DiagOperator(Hxx), Hilbert)
    return (; Hilbert, H)
end

# Define parameters
J = 1.0       # Interaction strength
h = 1.0       # Transverse field strength
lattice_size = 4  # Number of spins
periodic = false  # No periodic boundary conditions
# Define the Hamiltonian
(;Hilbert,H) = transverse_field_ising(lattice_size, h, J; periodic)
```
The Hamiltonian is split into a diagonal and an offdiagonal part. The diagonal $H_{xx}$ can be an arbitray function of the configuration $x$. The offdiagonal is given by the `moves` and `offdiagElements` arrays. The `moves` array contains the indices of the spins that are flipped, while the `offdiagElements` array contains the corresponding weights for each move.

!!! note
    The `moves` can also be given by Integers, i.e. `Int8[0,0,...,-1,1,0,...,0]` for a term $\sigma_i^- \sigma_{i+1}^+$. This will be slower, but allows for more complex moves.

## Running the Simulation

We now proceed to solve the Hamiltonian. It is instructive to consider a single walker first. As a guiding wavefunction, we use the simplest one, an equal weight superposition of all configurations $\psi(x) =1$.
!!! note
    Observers store simulation data. You can also initialize an Observer by passing a file name, i.e. `ConfigObserver(filename,startConfig, NSteps, NWalkers)` which will store everything in the file in HDF5 format by using memory mapping. Note that if this file is deleted when the Observer is used (i.e. in the simulation), julia will crash.

```@example TFI
NSteps = 500  # Number of Monte Carlo steps
NStepsEquil = 30  # Number of Monte Carlo steps for equilibration
NWalkers = 1  # Number of walkers
dtau = 0.1    # imaginary time propagation step

startConfig = BosonConfig(Hilbert) # creates an initial configuration where all occupation numbers are 0

problem = GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ = EqualWeightSuperposition(), H, Hilbert)
Observer = ConfigObserver(startConfig, NSteps, NWalkers) # Observer to measure the energy and configurations 
runGFMC!(problem, NoObserver(), NStepsEquil) #run for NStepsEquil steps without observing to equilibrate
runGFMC!(problem, Observer, NSteps) #run for NSteps steps
```
## Evaluating the Results
Let's see the results. For open boundary conditions at the critical point $h=J$, we may compare to the exact solution. 
```@example TFI
using CairoMakie, Statistics, Random
Random.seed!(123) # for reproducibility
function E_critPoint_exact(L, h=1, periodic=false)
    (!periodic && h==1) || return NaN 
    
    return 1 - csc(pi / (2 * (2 * L + 1)))
end

MaxProjection = 40
energies = getEnergies(Observer, MaxProjection) 
tau = (0:MaxProjection-1) * dtau

let 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$τ$ (imaginary time)", ylabel=L"Energy$$")

    lines!(ax, tau, energies, label=L"$\psi(x)=1$")

    hlines!(ax, E_critPoint_exact(lattice_size,h,periodic), color=:black, label=L"Exact$$", linestyle=:dash)

    axislegend(ax, position=:rt)
    fig
end
```
Run the GFMC simulation a few times (with different seeds) and see what happens:
- At $\tau=0$ (for a single walker!), we will always obtain the variational energy of the guiding wavefunction.
- The energy initially decreases with the number of projections.
- However, for larger $\tau$, the energy is strongly affected by statistical fluctuations. It may either increase or go below the exact energy.
- We may quantify this statistical error by looking at the standard deviation of the energy upon performing several runs (see below).
- Provided a large enough $\tau$, such that the imaginary time projection is converged, the energy will be within the statistical error of the exact energy.

## Using Variational Wavefunctions

The issue of errorbars growing exponentially with $\tau$ can be significantly alleviated by using more Walkers, i.e. setting `NWalkers = 10` in the example above, or by using a more sophisticated guiding wavefunction. For example, it is easy to implement a short ranged Jastrow wavefunction which correlates nearest neighbors spins.
!!! warning
    You must always implement the logarithm of the wavefunction, i.e. $\log \psi(x)$.

```@example TFI
function ShortRangeJastrow(x)
    res = 0.
    for i in eachindex(x)[1:end-1]
        res += 0.5x[i]*x[i+1]
    end
    return res
end
logψ = NaiveFunction(ShortRangeJastrow)
NWalkers = 10
P = ProblemEnsemble([GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ, H, Hilbert) for i in 1:10])
Observers_jastrow = [ConfigObserver(startConfig, NSteps, NWalkers) for _ in 1:10]

runGFMC!(P, NoObserver(),200,logger=nothing) #equilibrate
runGFMC!(P, Observers_jastrow,NSteps,logger=nothing)

energies_jastrow = [getEnergies(Observer, MaxProjection) for Observer in Observers_jastrow]

let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$τ$ (imaginary time)", ylabel=L"Energy$$")

    lines!(ax, tau, energies, label=L"$\psi(x)=1,\ N_w=1$")

    lines!(ax, tau, mean(energies_jastrow), label=L"shortrange Jastrow $N_w=%$NWalkers$")
    band!(ax, tau, mean(energies_jastrow) - std(energies_jastrow), mean(energies_jastrow) + std(energies_jastrow), color=(:green, 0.2))

    hlines!(ax, E_critPoint_exact(lattice_size), color=:black, label=L"Exact$$", linestyle=:dash)

    axislegend(ax, position=:rt)
    fig
end
```
!!! tip
    There is already an efficient implementation of the Jastrow function. Simply use it as:
    ```
    vij_jastrow = zeros(Float32,lattice_size,lattice_size)
    for i in axes(vij_jastrow, 1)[1:end-1]
        vij_jastrow[i,i+1]=vij_jastrow[i+1,i] = 0.5
    end

    logψJastrow = Jastrow(
        zeros(Float32,lattice_size),
        vij_jastrow,
    )
    ```
    It is advised to use variational Monte Carlo to optimize a variational wavefunction before using it in GFMC. A particularly useful package is [Netket](https://netket.readthedocs.io/en/stable/), which may be called from Julia via [PyCall](https://github.com/JuliaPy/PyCall.jl).


## Observables diagonal in the computational basis
Observables which are diagonal in the computational basis (i.e. they can be expressed as a function of the occupation numbers) are simple to determine. Using the `ConfigObserver`, which records the configuration of the walkers at each step, we can compute them cheaply after the simulation. 
To do this, we only need to use the function `getObs_diagonal`. Let's consider the average occupation number as an example.
```@example TFI
mProj = 40  # Number of projection steps
Observable = OccupationNumber(lattice_size)

n_avg = [stack(getObs_diagonal(O,Observable,1:mProj)) for O in Observers_jastrow]

n_avg_mean = mean(n_avg)
n_avg_std = std(n_avg)

let 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$τ$ (imaginary time)", ylabel=L"$\langle S^z_i\rangle$")
    tau =( 0:mProj-1) * dtau
    for i in 1:lattice_size
        lines!(ax, tau, n_avg_mean[i,:])
        band!(ax, tau, n_avg_mean[i,:] - n_avg_std[i,:], n_avg_mean[i,:] + n_avg_std[i,:], alpha = 0.5)
    end
    fig
end
```

Here, `OccupationNumber` is a pre-defined observable. Correlation functions such as $\langle S^z_i S^z_j \rangle$ are common enough to also come pre-defined: [`SpinCorrelations`](@ref) is a `@turbo`-accelerated implementation that only computes the upper triangle of the correlation matrix (use [`get_matrix_from_tri`](@ref) to expand its output back into a full matrix).

It is also relatively simple to implement your own observables by defining a new subtype `O` of [`AbstractObservable`](@ref): you need `obs(::O)` returning a preallocated output buffer, a callable `(::O)(out, config)` that fills `out` given a configuration, and `Base.copy(::O)`. The `OccupationNumber` and `SpinCorrelations` source files are good templates to start from.

!!! note "Two ways to measure diagonal observables"
    So far we used [`ConfigObserver`](@ref) together with `getObs_diagonal` to compute observables *after* the simulation, from the stored configurations. This is the most flexible option: any observable can be evaluated at any projection step, purely in post-processing. Its cost is memory/storage — every walker configuration at every observation step must be kept.

    For long-running simulations, it is usually preferable to accumulate observables *during* the run instead, using a [`LazyObservableAccumulator`](@ref). This fixes the set of observables and projection steps ahead of time, but avoids storing full configuration histories, which lets you scale to far more walkers/steps. This is demonstrated next.

## Efficient observable accumulation for large runs

For production runs, [`LazyObservableAccumulator`](@ref) is the recommended way to measure diagonal observables such as [`OccupationNumber`](@ref) and [`SpinCorrelations`](@ref): it projects each observable to the requested imaginary times on the fly and only stores the resulting (small) numerator/denominator arrays, rather than every walker configuration.

We first run a few epochs without observing to get a good estimate of the average walker weight — this improves numerical stability of the accumulation and can replace an initial equilibration run. We also use [`BatchMultiThreaded`](@ref), a Polyester-based parallelization scheme that tends to outperform the default `MultiThreaded` scheme when many CPU cores are available (see the "Parallelization" section in [Advanced Usage](Advanced.md) for details).

```julia
NWalkers = 10
parallelization = BatchMultiThreaded(Threads.nthreads())
problem = GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ, H, Hilbert, parallelization)

mean_TotalWeights, w_avg_estimate = GFMC.estimate_weights_continuousTime!(problem; Nepochs=4, Nsamples=1000, verbose=true, logger=nothing)
```

The set of accumulators is combined into a single `CombinedObserver`, so that each of them is updated at every observation step. `BasicAccumulator` provides the buffers shared by all `LazyObservableAccumulator`s, and stores the energy itself.

```julia
outfile = "observables.h5"
NSteps_total = 2000  # total number of simulation steps
num_bins = 32         # number of bins used for the accumulation (see "Advanced Usage" for error estimation via binning)
projection_values = 0:5:mProj-1 # imaginary-time projections at which the observables are stored; storing fewer steps is cheaper

BAcc = BasicAccumulator(outfile, mProj, NWalkers; weight_normalization = mean_TotalWeights, num_bins, bin_elements = NSteps_total ÷ num_bins + 1)
OccNum_Acc = LazyObservableAccumulator(outfile, BosonConfig(Hilbert), OccupationNumber(lattice_size), BAcc, projection_values, NWalkers, Threads.nthreads())
SpinCorr_Acc = LazyObservableAccumulator(outfile, BosonConfig(Hilbert), SpinCorrelations(lattice_size), BAcc, projection_values, NWalkers, Threads.nthreads())
Observer = CombinedObserver((BAcc, OccNum_Acc, SpinCorr_Acc))

CT = ContinuousTimePropagator(dtau, w_avg_estimate) # use the estimated average weight to reduce numerical error
runGFMC!(problem, Observer, NSteps_total; Propagator = CT, parallelization)
```

After the run, the numerators and denominators for each observable are available both on the accumulators themselves and in `outfile` (HDF5). See [Advanced Usage](Advanced.md) for how to bin the data, estimate error bars, and read it back from disk, as well as for the buffered [`ObservableAccumulator`](@ref) — a middle ground that keeps observable buffers in memory across a few projection steps, at extra memory cost, to avoid recomputing expensive observables.
