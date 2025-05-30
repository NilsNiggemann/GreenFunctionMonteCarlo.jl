include("utils.jl") 

function E_critPoint_exact(L)
    return 1 - csc(pi / (2 * (2 * L + 1)))
end
σz(n::Bool) = (1 - 2 * n)
σz(i, conf::AbstractArray) = σz(conf[i])

using GreenFunctionMonteCarlo.LinearAlgebra
NSites = 2
NSteps = 1500
mProj = 20
NWalkers = 2

RNG = StableRNG(1234)

Hilbert = BosonHilbertSpace(NSites, HardCoreConstraint())
moves = eachcol(Bool.(I(NSites))) # each move flips a single spin
offdiagElements = -ones(NSites)
Hxx = DiagOperator(x-> sum(σz(i, x) * σz(i + 1, x) for i in eachindex(x)[1:end-1]))

H = localOperator(moves, offdiagElements, Hxx, Hilbert)

config = BosonConfig(Hilbert)
rand!(RNG,config)
logψ = GFMC.EqualWeightSuperposition()
CT = ContinuousTimePropagator(0.1,w_avg_estimate  = -2.09)

prob = GFMCProblem(config, NWalkers, CT; logψ, H, Hilbert)

outfile = tempname()

BasicAccumulatorFile = GFMC.BasicAccumulator(outfile, mProj, NWalkers)
ObsAccumulatorFile = GFMC.ObservableAccumulator(outfile,OccupationNumber(zeros(NSites)),BasicAccumulatorFile, mProj, NWalkers, Threads.nthreads())

BObs = GFMC.BasicObserver(outfile, NSteps, NWalkers)
CObs = GFMC.ConfigurationObserver(outfile, config, NSteps, NWalkers)
Observer = GFMC.CombinedObserver((BObs, CObs,BasicAccumulatorFile, ObsAccumulatorFile))
runGFMC!(prob, NoObserver(), 200; rng = RNG)
runGFMC!(prob, Observer, NSteps; rng = RNG)
##
Gnps = GFMC.precomputeNormalizedAccWeight(BObs.TotalWeights,2mProj)
n_avg = stack(getObs_diagonal(Gnps,CObs.SaveConfigs,BObs.reconfigurationTable,OccupationNumber(zeros(NSites)),1:mProj))

n_avg_direct = ObsAccumulatorFile.Obs_numerator ./ ObsAccumulatorFile.Obs_denominator'

