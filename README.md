# GreenFunctionMonteCarlo.jl
[![build](https://github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl/workflows/CI/badge.svg)](https://github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl/actions?query=workflow%3ACI) [![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://NilsNiggemann.github.io/GreenFunctionMonteCarlo.jl/dev)
<!-- Tidyverse lifecycle badges, see https://www.tidyverse.org/lifecycle/ Uncomment or delete as needed. -->
![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->

<!-- [![GitHub commits since tagged version](https://img.shields.io/github/commits-since/NilsNiggemann/GreenFunctionMonteCarlo.jl/v0.0.1.svg?style=social&logo=github)](https://github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl) -->
<!-- [![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://NilsNiggemann.github.io/GreenFunctionMonteCarlo.jl/stable) -->

<!-- Aqua badge, see test/runtests.jl -->
<!-- [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) -->

<!-- travis-ci.com badge, uncomment or delete as needed, depending on whether you are using that service. -->
<!-- [![Build Status](https://travis-ci.com/NilsNiggemann/GreenFunctionMonteCarlo.jl.svg?branch=master)](https://travis-ci.com/NilsNiggemann/GreenFunctionMonteCarlo.jl) -->
<!-- NOTE: Codecov.io badge now depends on the token, copy from their site after setting up -->
<!-- Documentation -- uncomment or delete as needed -->
<!-- START README (DO NOT DELETE THIS LINE!) -->
## Overview

`GreenFunctionMonteCarlo.jl` is a Julia package designed for performing Green Function Monte Carlo (GFMC) simulations on lattice models.

Presently, this package treats only Hamiltonians that are free of the sign problem, i.e. whose elements satisfy
```math

H_{x, x'} \leq 0 \quad \forall x \neq x'
```
where $H_{x, x'}$ is the matrix element of the Hamiltonian between two configurations (spins or bosons) $x$ and $x'$. 

## Features
- [x] Efficient implementation of the fixed-walker formalism for general systems of interacting bosons [\[1\]](#references).
- [x] Flexible framework for implementation of other algorithms.
- [x] Constraints, such as hardcore boson constraint.
- [x] Interface for the efficient implementation of arbitrary constraints.
- [ ] Measure arbitrary observables in the diagonal basis.
- [ ] Allow more direct implementation of spin models.
- [ ] Allow for easy conversions of Hamiltonians and variational wavefunctions from netket.

## Limitations
- The package is currently in the experimental stage and exported functionality may change or break in the future.
- Currently, there is no implementation of the fixed-node approximation [\[2\]](#references) for Hamiltonians with the sign problem.

## Installation

To install the package, use the Julia package manager:

```julia
using Pkg
Pkg.add(url = "https://github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl.git")
```

## Documentation

For detailed documentation, including API references and advanced usage, visit the [documentation site](https://NilsNiggemann.github.io/GreenFunctionMonteCarlo.jl/dev).

## Contributing

Contributions are welcome! If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request on the [GitHub repository](https://github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl).

## License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl/blob/master/LICENSE) file for details.


## References:
- [1] [Buonaura, M. & Sorella, S. Green's function Monte Carlo method for lattice fermions. Phys. Rev. B 57, 11446 (1998).](https://doi.org/10.1103/PhysRevB.57.11446)
- [2] [Becca, F. & Sorella, S. *Quantum Monte Carlo Approaches for Correlated Systems*. Cambridge University Press; 2017.](https://doi.org/10.1017/9781316417041)