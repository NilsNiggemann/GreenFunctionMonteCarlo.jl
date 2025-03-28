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
The Hamiltonian is split into a diagonal and an offdiagonal part. The diagonal $H_{xx\prime}$ can be an arbitray function of the configuration $x$. The offdiagonal is given by the `moves` and `offdiagElements` arrays. The `moves` array contains the indices of the spins that are flipped, while the `offdiagElements` array contains the corresponding weights for each move.

!!! note
    The `moves` can also be given by Integers, i.e. `Int8[0,0,...,-1,1,0,...,0]` for a term $\sigma_i^- \sigma_{i+1}^+$. This will be slower, but allows for more complex moves.

## Running the Simulation

We now proceed to solve the Hamiltonian. It is instructive to consider a single walker first. As a guiding wavefunction, we use the simplest one, an equal weight superposition of all configurations $\psi(x) =1$.

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
tau = 0:MaxProjection-1 * dtau

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