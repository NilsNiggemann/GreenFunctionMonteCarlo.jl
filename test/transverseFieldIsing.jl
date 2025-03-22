using GreenFunctionMonteCarlo
import GreenFunctionMonteCarlo as GFMC
##
using Random
using StableRNGs
import GreenFunctionMonteCarlo.SmallCollections as SC
# https://dmrg101-tutorial.readthedocs.io/en/latest/tfim.html
function getTFI(Nsites, h, J, rng=Random.default_rng())
    Hilbert = BosonHilbertSpace(Nsites, HardCoreConstraint())
    moves = [Bool[0 for _ in 1:Nsites] for _ in 1:Nsites]
    weights = zeros(Float64, Nsites)

    for i in 1:Nsites
        moves[i][i] = true
        weights[i] = -h
    end

    function σz(i, conf)
        return 0.5*(1-2*conf[i])
    end
    function Hxx(conf)
        
        E = -J * sum(σz(i,conf)*σz(i+1,conf) for i in 1:Nsites-1)
        return E
    end
    H = localOperator(moves, weights, GFMC.DiagOperator(Hxx), Hilbert)
    return (; Hilbert, H)
end
function E_critPoint(L)
    return 1 - csc(pi / (2 * (2 * L + 1)))
end

RNG = StableRNG(1234)
NSites = 20
(Hilbert,H) = getTFI(NSites,1,1)

config = BosonConfig(Hilbert)

rand!(RNG,config)

# logψ = Jastrow(config,Float64)
# params = get_params(logψ)
# rand!(RNG,params)
# params .*= 1e-3
# logψ.v_ij .= GFMC.LinearAlgebra.Symmetric(logψ.v_ij)
logψ = EqualWeightSuperposition()

NWalkers = 10000
NSteps = 1000
CT = ContinuousTimePropagator(0.1,-12.0)
# logger = GFMC.SimpleLogger(10)
logger = GFMC.ProgressBarLogger(dt=0.01)
prob = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert,parallelization = GFMC.SingleThreaded())
ConfSaver = ConfigObserver(config, NSteps,NWalkers)

runGFMC!(prob, NoObserver(), 100;rng= RNG)

runGFMC!(prob, ConfSaver, NSteps;rng= RNG, logger=ProgressBarLogger())
println("E_ex = ", E_critPoint(NSites))
(;SaveConfigs,TotalWeights,energies,reconfigurationTable) = ConfSaver