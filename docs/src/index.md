# GreenFunctionMonteCarlo
## Overview

`GreenFunctionMonteCarlo.jl` is a Julia package designed for performing Green Function Monte Carlo (GFMC) simulations on lattice models.

Presently, this package treats only Hamiltonians that are free of the sign problem, i.e. whose elements satisfy
```math

H_{x, x'} \leq 0 \quad \forall x \neq x'
```
where $H_{x, x'}$ is the matrix element of the Hamiltonian between two configurations (spins or bosons) $x$ and $x'$. 

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