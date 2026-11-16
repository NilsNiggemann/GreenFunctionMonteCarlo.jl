```@meta
CurrentModule = GreenFunctionMonteCarlo
```

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
    There is already an efficient implementation of the Jastrow function, used in the [Accumulating Observables On-The-Fly](@ref) section below. It is advised to use variational Monte Carlo to optimize a variational wavefunction before using it in GFMC. A particularly useful package is [Netket](https://netket.readthedocs.io/en/stable/), which may be called from Julia via [PyCall](https://github.com/JuliaPy/PyCall.jl).


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

mean_TotalWeights, w_avg_estimate = estimate_weights_continuousTime!(problem; Nepochs=4, Nsamples=1000, verbose=true, logger=nothing)
```

The set of accumulators is combined into a single `CombinedObserver`, so that each of them is updated at every observation step. `BasicAccumulator` provides the buffers shared by all `LazyObservableAccumulator`s, and stores the energy itself.

!!! tip
    While the approach above is very convenient, for big simulations, it may not be feasible to store all configurations as the output file may become too large. For this case, it is also possible to use accumulators, such as [`BasicAccumulator`](@ref) and [`LazyObservableAccumulator`](@ref). Accumulators compute the imaginary time projection of the observable at every step of the simulation, thereby saving a lot of storage. `BasicAccumulator` contains all the essential information to allow for projection during the run, while `LazyObservableAccumulator` may be used to compute observables. 
    To combine several accumulators, you can use `CombinedObserver`. The following section shows this in practice.

## Accumulating Observables On-The-Fly

For longer simulations, storing the full configuration at every step (as `ConfigObserver` does) becomes wasteful once you only care about a handful of expectation values. Instead, [`BasicAccumulator`](@ref) and [`LazyObservableAccumulator`](@ref) accumulate the imaginary-time projection of each observable step by step and write directly to an HDF5 file, so the memory and storage footprint no longer grows with the number of steps. This comes at the cost of having to fix the maximum projection depth `m_proj` before running the simulation. See [Advanced Usage](Advanced.md) for a more detailed discussion, including numerical stability caveats and the buffered [`ObservableAccumulator`](@ref).

We reuse the `Hilbert`, `H` and `logψ` (short-range Jastrow) already defined above, but this time with `NWalkers = 10` and a fixed number of bins, so that we can later estimate error bars.
```@example TFI
using GreenFunctionMonteCarlo.HDF5

function adjust_NSteps_for_bins(NSteps, num_bins)
    bin_elements = NSteps ÷ num_bins + 1
    return num_bins * bin_elements - 1
end

NWalkers      = 10
m_proj        = 40   # maximal projection order used by BasicAccumulator/LazyObservableAccumulator
obs_projection_stepsize = 2 # skip every other projection step for the (cheaper) observables
num_bins      = 32   # number of bins used for the accumulation, e.g. for error estimation via binning
NSteps        = adjust_NSteps_for_bins(500, num_bins)
outfile       = tempname() * ".h5" # accumulators write their running sums directly to this file
nothing # hide
```

The `BasicAccumulator` stores the running numerator/denominator of the energy estimator, split into `num_bins` bins. `LazyObservableAccumulator` reuses this bookkeeping to accumulate further observables; here we measure the occupation numbers and the spin-spin correlations. Both accumulators are packed into a `CombinedObserver`, which updates each of them at every observation step; `BasicAcc` must come first, since the other accumulators read its internal buffers.
```@example TFI
problem = GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ, H, Hilbert)
mean_TotalWeights, w_avg_estimate = estimate_weights_continuousTime!(problem; Nepochs=4, Nsamples=1000, verbose=true, logger=nothing)
problem = GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau, w_avg_estimate); logψ, H, Hilbert)

BasicAcc = BasicAccumulator(outfile, m_proj, NWalkers; num_bins, bin_elements = NSteps ÷ num_bins + 1, weight_normalization = mean_TotalWeights)

OccNum = OccupationNumber(lattice_size)
SpinCorr = SpinCorrelations(lattice_size)
projection_values = 0:obs_projection_stepsize:(m_proj-1)
OccNumAcc = LazyObservableAccumulator(outfile, startConfig, OccNum, BasicAcc, projection_values, NWalkers, 1)
SSAcc = LazyObservableAccumulator(outfile, startConfig, SpinCorr, BasicAcc, projection_values, NWalkers, 1)

Observer = CombinedObserver((BasicAcc, OccNumAcc, SSAcc))

runGFMC!(problem, Observer, NSteps);
nothing # hide
```

### Reading the results back from the HDF5 file

Once the simulation is done, the accumulated numerators and denominators live in `outfile`. As described in [Advanced Usage](Advanced.md), the recommended way to post-process them is to read them back from the file (rather than reaching into the live accumulator objects) and bin them together. [`get_energy_from_accumulator_bunching`](@ref) and [`get_obs_from_accumulator_bunching`](@ref) can read directly from the HDF5 file, given its path and (for observables) the name they were stored under:
```@example TFI
function readObs(file, bunching)
    Energy = get_energy_from_accumulator_bunching(file, bunching)
    OccupationNumbers = get_obs_from_accumulator_bunching(file, "OccupationNumber", bunching)
    SzSz = get_obs_from_accumulator_bunching(file, "SpinCorrelations", bunching)

    m_values = h5read(file, "OccupationNumber_m_values")
    return (; Energy, OccupationNumbers, SzSz, m_values)
end

bunching = num_bins ÷ 8 # bin the num_bins accumulated bins down to a handful of groups for error bars
Obs = readObs(outfile, bunching)
nothing # hide
```
`Obs.Energy`, `Obs.OccupationNumbers` and `Obs.SzSz` are each vectors with one entry per (bunched) bin, which we can average over to get the expectation values together with an error estimate.
```@example TFI
tau_energy = (0:m_proj-1) .* dtau
energy_mean = mean(Obs.Energy) .* lattice_size
energy_std = std(Obs.Energy) .* lattice_size

let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"$τ$ (imaginary time)", ylabel = L"Energy$$")

    lines!(ax, tau_energy, energy_mean, label = L"GFMC$$")
    band!(ax, tau_energy, energy_mean .- energy_std, energy_mean .+ energy_std, color = (:blue, 0.2))
    hlines!(ax, E_critPoint_exact(lattice_size, h, periodic), color = :black, linestyle = :dash, label = L"Exact$$")

    axislegend(ax, position = :rt)
    fig
end
```
```@example TFI
n_mean = mean(Obs.OccupationNumbers)
n_std = std(Obs.OccupationNumbers)
tau_obs = collect(projection_values) .* dtau

let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"$τ$ (imaginary time)", ylabel = L"$\langle n_j \rangle$")
    for j in 1:lattice_size
        lines!(ax, tau_obs, n_mean[j, :], label = L"j=%$j")
        band!(ax, tau_obs, n_mean[j, :] .- n_std[j, :], n_mean[j, :] .+ n_std[j, :], alpha = 0.2)
    end
    axislegend(ax, position = :rt)
    fig
end
```
```@example TFI
S1j_inds = findall(==(1), SpinCorr.i_inds)   # rows for pairs (1,j), ordered j = 1..lattice_size
SzSz_mean = mean(Obs.SzSz)[S1j_inds, :]
SzSz_std = std(Obs.SzSz)[S1j_inds, :]

let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"$τ$ (imaginary time)", ylabel = L"$\langle S_1^z S_j^z \rangle$")
    for j in 1:lattice_size
        lines!(ax, tau_obs, SzSz_mean[j, :], label = L"j=%$j")
        band!(ax, tau_obs, SzSz_mean[j, :] .- SzSz_std[j, :], SzSz_mean[j, :] .+ SzSz_std[j, :], alpha = 0.2)
    end
    axislegend(ax, position = :rt)
    fig
end
```

!!! note
    Here we only bin the `num_bins` accumulated bins down once, with a single `bunching` value, to get *some* error estimate. For a real production run, you should instead check that this error estimate is stable as the bunching size grows (i.e. that the bins are sufficiently decorrelated) before trusting it — see the warning in the "Example evaluation" section of [Advanced Usage](Advanced.md).

## Observables off-diagonal in the computational basis

Observables such as the transverse magnetization $\sum_i \sigma^x_i$ change the configuration when applied, so they cannot be read off from the walker configurations, and `getObs_diagonal` does not apply. Instead, they are measured by **forward walking**: every `cadence` steps, the operator is applied once to a snapshot of the running walker population. The resulting sub-population is then propagated for `mProj` further steps with the ordinary dynamics, and is normalized against the main population's own weight growth. For depth $p=0$ this gives the biased *mixed estimator* $\langle \psi_T | O | \psi_0 \rangle$; as $p$ grows, the estimate converges to the true ground-state expectation value $\langle O \rangle$, in the same way that the energy converges with the projection order.

Three ingredients are needed:

1. An off-diagonal observable, built with [`offdiagonalObservable`](@ref) from a list of `moves` (same format as for the Hamiltonian) and a weight function `weightfunc(idx, ψratio, x)`. Here, `ψratio` is the ratio $\psi(x')/\psi(x)$ of the guiding function for move `idx`, and the weight must be non-negative. For $\sigma^x_i$ the matrix element is $1$, so the weight is just `ψratio`. If a move can be forbidden (e.g. by a hardcore constraint), return `0` for it.
2. A [`ForwardWalkingAccumulator`](@ref), which needs the same `logψ`, `H`, `Hilbert` and propagator as the main simulation, the forward depth `mProj`, the `cadence` and the sub-population size `Nsub` (which must equal `NWalkers`).
3. The same [`BasicAccumulator`](@ref) that is passed to the `CombinedObserver` must also be given to the `ForwardWalkingAccumulator`, and it must come *first* in the `CombinedObserver`. It supplies the reference-population weights used for normalization, and the bin bookkeeping. It requires `mProj < 2 * m_proj`.

!!! note
    `offdiagonalObservable` and `ForwardWalkingAccumulator` are not exported, so we import them explicitly.

We reuse the setup of the previous section (`H`, `Hilbert`, `startConfig`, `m_proj`, `num_bins`, `NSteps`, `NWalkers`). The forward-walking sub-populations need a guiding function with an efficient update buffer, so instead of the `NaiveFunction` used above we use the built-in [`Jastrow`](@ref) with the same nearest-neighbor correlations.
```@example TFI
using GreenFunctionMonteCarlo: offdiagonalObservable, ForwardWalkingAccumulator,
    get_observable_from_accumulator_bunching

# σˣ_i flips site i, with matrix element 1, so the weight is the bare guiding-function ratio.
# Since all moves are included, the result is the sum Σ_i ⟨σˣ_i⟩ over all sites.
sx_moves = [Bool[i == j for j in 1:lattice_size] for i in 1:lattice_size]
SigmaX = offdiagonalObservable(sx_moves, (idx, ψratio, x) -> ψratio, Hilbert)
```
```@example TFI
mProj_fw = 20     # number of forward-walking steps; the estimate is obtained for p = 0, ..., mProj_fw
cadence = 2       # seed a new sub-population every 2 steps
outfile_fw = tempname() * ".h5"

v_ij = zeros(lattice_size, lattice_size) # Jastrow couplings, nonzero for nearest neighbors only
for i in 1:lattice_size-1
    v_ij[i, i+1] = v_ij[i+1, i] = 0.5
end
logψ_fw = Jastrow(zeros(lattice_size), v_ij)

problem = GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ = logψ_fw, H, Hilbert)
mean_TotalWeights_fw, w_avg_fw = estimate_weights_continuousTime!(problem; Nepochs=4, Nsamples=1000, verbose=false, logger=nothing)
propagator = ContinuousTimePropagator(dtau, w_avg_fw)
problem = GFMCProblem(startConfig, NWalkers, propagator; logψ = logψ_fw, H, Hilbert)

BasicAcc_fw = BasicAccumulator(outfile_fw, m_proj, NWalkers; num_bins, bin_elements = NSteps ÷ num_bins + 1, weight_normalization = mean_TotalWeights_fw)
FWAcc = ForwardWalkingAccumulator(outfile_fw, SigmaX, BasicAcc_fw, startConfig, logψ_fw, H, Hilbert,
                                  propagator, mProj_fw, cadence, NWalkers)

runGFMC!(problem, CombinedObserver((BasicAcc_fw, FWAcc)), NSteps; logger=nothing);
nothing # hide
```
As for the diagonal accumulators, we bin the result to obtain error bars. The output has one vector of length `mProj_fw + 1` (the forward depth $p$) per bunched bin.
```@example TFI
sx_bins = get_observable_from_accumulator_bunching(FWAcc, bunching)
sx_mean = mean(sx_bins)
sx_std = std(sx_bins)
nothing # hide
```
For comparison, we diagonalize the $2^L$-dimensional Hamiltonian exactly (feasible for our small system), build the same operator $\sum_i \sigma^x_i$, and plot the estimate as a function of the forward depth.
```@example TFI
using LinearAlgebra

function exact_sigmax(L, h, J; periodic = false)
    D = 2^L
    Hmat = zeros(D, D)
    Ox = zeros(D, D)
    for s in 0:D-1
        spins = [1 - 2 * ((s >> (i - 1)) & 1) for i in 1:L] # σᶻ_i = ±1
        bonds = periodic ? L : L - 1
        Hmat[s+1, s+1] = -J * sum(spins[i] * spins[mod1(i + 1, L)] for i in 1:bonds)
        for i in 1:L
            t = xor(s, 1 << (i - 1)) # flip spin i
            Hmat[t+1, s+1] += -h
            Ox[t+1, s+1] += 1
        end
    end
    ψ0 = eigen(Symmetric(Hmat)).vectors[:, 1]
    return ψ0' * Ox * ψ0
end
sx_exact = exact_sigmax(lattice_size, h, J; periodic)

let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"forward depth $p$", ylabel = L"$\sum_i \langle \sigma^x_i \rangle$")
    p = 0:mProj_fw
    lines!(ax, p, sx_mean, label = L"GFMC$$")
    band!(ax, p, sx_mean .- sx_std, sx_mean .+ sx_std, color = (:blue, 0.2))
    hlines!(ax, sx_exact, color = :black, linestyle = :dash, label = L"Exact$$")
    axislegend(ax, position = :rb)
    fig
end
```
At $p=0$ we see the biased mixed estimator, which depends on the guiding wavefunction. With increasing $p$, the estimate moves towards the exact value and is compatible with it within the error bars (which grow with $p$, as for the energy). Convergence should be checked for convergence with `mProj_fw` in the same way as the projection order of the energy.

!!! tip
    - To measure a *momentum-resolved* observable, such as a structure factor $B(q)=\sum_j e^{iqj}\sigma^x_j$, the weights must stay non-negative for sampling. A `cos²` phase weighting is a convenient choice, see [`offdiagonalObservable`](@ref).
    - By default, a multi-move observable yields the *sum* over all its moves. Pass `resolve_moves = true` to `ForwardWalkingAccumulator` to obtain one estimate per move (an `(n_moves, mProj+1)` array per bin) instead, at `n_moves` times the cost.
    - Several accumulators can share one `BasicAccumulator` in a `CombinedObserver`; their sub-populations are then advanced together. Use a multithreaded `parallelization` to speed up seeding and propagation.
    - Forward walking costs roughly `fld(mProj, cadence) + 1` additional small populations per step, but memory does not grow with `NSteps`.
