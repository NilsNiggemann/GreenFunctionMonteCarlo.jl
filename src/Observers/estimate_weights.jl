"""
    estimate_weights_continuousTime!(prob; Nepochs=5, Nsamples=100, mProj=50)

Estimate weights in a continuous-time Monte Carlo simulation. Useful for reducing the floating point errors for accumulators. The results may be passed to 
`ContinuousTimePropagator(tau,w_avg_estimate)` and `BasicAccumulator(args...;weight_normalization)`

# Arguments
- `prob`: The problem instance or data structure containing the simulation setup.
- `Nepochs`: (Optional) Number of epochs to run the estimation. Default is 5.
- `Nsamples`: (Optional) Number of samples per epoch. Default is 100.
- `mProj`: (Optional) Number of projection steps or iterations. Default is 50.

# Description
This function performs an in-place estimation of weights for a given problem using a continuous-time Monte Carlo approach. The process iterates over a specified number of epochs, drawing samples and projecting as configured by the keyword arguments.

# Returns
- Modifies `prob` in place which helps to equilibrate
- Returns a rough estimate of the mean total weight and CT.w_avg_estimate.
"""
function estimate_weights_continuousTime!(prob;Nepochs=5,Nsamples=100,mProj = 50,verbose = false)
    Observer = BasicObserver(Nsamples, NWalkers(prob.Walkers))
    CT = prob.Propagator
    mean_TotalWeights = 1.
    for epoch in 1:Nepochs
        runGFMC!(prob, Observer, Nsamples; Propagator = CT)
        mean_TotalWeights = Statistics.mean(Observer.TotalWeights)
        E0tau = getEnergies(Observer.TotalWeights, Observer.energies, mProj)
        E0_ind = findfirst(<=(0), diff(E0tau))
        isnothing(E0_ind) && (E0_ind = mProj)
        E0 = E0tau[E0_ind]
        CT = ContinuousTimePropagator(CT.dτ,-E0)
        if verbose 
            @info "" epoch mean_TotalWeights E0
        end
    end
    return mean_TotalWeights, CT.w_avg_estimate
end