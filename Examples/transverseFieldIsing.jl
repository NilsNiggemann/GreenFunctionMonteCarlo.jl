using GreenFunctionMonteCarlo
import GreenFunctionMonteCarlo as GFMC
using CairoMakie, MakieHelpers
using Statistics
##
using Random
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
        return (1-2*conf[i])
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

# logψ = Jastrow(config,Float64)
# params = get_params(logψ)
# rand!(RNG,params)
# params .*= 1e-3
# logψ.v_ij .= GFMC.LinearAlgebra.Symmetric(logψ.v_ij)

function makeRuns(prob,ConfSaver,Nruns)

    res = Vector{typeof(ConfSaver)}(undef,Nruns)
    Threads.@threads for i in 1:Nruns
        newProb = deepcopy(prob)
        rand!.(newProb.WE.Configs)

        newConfSaver = deepcopy(ConfSaver)
        runGFMC!(newProb, NoObserver(), 100,logger = NoLogger())
        runGFMC!(newProb, newConfSaver, NSteps,logger = NoLogger())
        res[i] = newConfSaver
    end
    return res
end


NSites = 20
(Hilbert,H) = getTFI(NSites,1,1)

config = BosonConfig(Hilbert)

logψ = EqualWeightSuperposition()

NWalkers = 28*20
NSteps = 8000
CT = ContinuousTimePropagator(0.1,-12.0)
prob = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert,parallelization = GFMC.SingleThreaded())
ConfSaver = ConfigObserver(config, NSteps,NWalkers)

results = makeRuns(prob,ConfSaver,20)
# (;SaveConfigs,TotalWeights,energies,reconfigurationTable) = ConfSaver
##
energies = getEnergies.(results,100)

let 
    fig = Figure()
    ax = Axis(fig[1, 1])
    errlines!(ax,mean(energies),std(energies),label="E")
    hlines!(ax,[E_critPoint(NSites),],label="E_ex")
    fig
    
end