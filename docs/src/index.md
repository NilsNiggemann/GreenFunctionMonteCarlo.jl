# GreenFunctionMonteCarlo
## Overview

`GreenFunctionMonteCarlo.jl` is a Julia package designed for performing Green Function Monte Carlo (GFMC) simulations on lattice models.

Presently, this package treats only Hamiltonians that are free of the sign problem, i.e. whose elements satisfy
```math

H_{x, x'} \leq 0 \quad \forall x \neq x'
```
where $H_{x, x'}$ is the matrix element of the Hamiltonian between two configurations (spins or bosons) $x$ and $x'$. 
## Installation

To install the package, use the Julia package manager:

```julia
using Pkg
Pkg.add(url = "https://github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl.git")
```

## Quick usage example: 
```julia
using GreenFunctionMonteCarlo, LinearAlgebra
NSites = 3
Nwalkers = 10
Hilbert = BosonHilbertSpace(NSites, HardCoreConstraint())
moves = Bool.(I(NSites)) # each move flips a single spin
offdiagElements = -ones(NSites)
H = localOperator(moves, offdiagElements, DiagOperator(x->0), Hilbert)

prob = GFMCProblem(BosonConfig(Hilbert), Nwalkers, ContinuousTimePropagator(0.1); logψ = EqualWeightSuperposition(), H, Hilbert)
Observer = ConfigObserver(startConfig, NSteps, NWalkers) # Observer to measure the energy and configurations
runGFMC!(problem, NoObserver(), NStepsEquil) #run for NStepsEquil steps without observing to equilibrate
runGFMC!(problem, Observer, NSteps) #run for NSteps steps
```

## Contributing

Contributions are welcome! If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request on the [GitHub repository](https://github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl).

## License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl/blob/master/LICENSE) file for details.


## References:
- [1] [Buonaura, M. & Sorella, S. Green's function Monte Carlo method for lattice fermions. Phys. Rev. B 57, 11446 (1998).](https://doi.org/10.1103/PhysRevB.57.11446)
- [2] [Becca, F. & Sorella, S. *Quantum Monte Carlo Approaches for Correlated Systems*. Cambridge University Press; 2017.](https://doi.org/10.1017/9781316417041)
