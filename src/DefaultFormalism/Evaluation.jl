function precomputeNormalizedAccWeight(weights,PMax;normalize=true)
    bn = weights
    meanweight = normalize ? Statistics.mean(bn) : 1
    
    Gnp = zeros(length(bn),PMax)

    for n in axes(Gnp,1)
        Gnp[n,1] = !iszero(bn[n]) #zero projection order
        # Gnp[n,1] = 1 #zero projection order
        PMax < 2 && continue
        Gnp[n,2] = bn[n]/meanweight # first projection order
    end
    for p in 3:PMax
        for n in p:length(bn)
            # Gnp[n,p] = Gnp[n,p-1]*bn[n-p]/meanweight
            Gnp[n,p] = Gnp[n-1,p-1]*Gnp[n,2]
        end
    end
    return Gnp
end

"""
    getEnergies(weights, localEnergies, PMax; kwargs...)

Compute the energies based on the provided weights and local energies, for the projection steps `p = 0,1,2,...,PMax`.

# Arguments
- `weights`: A collection of weights associated with the samples.
- `localEnergies`: A collection of local energy values corresponding to the samples.
- `PMax`: An integer indicating the maximum projection order.
# Keyword arguments
- `Gnp`: optionally provided precomputed normalized accumulated weights.

# Returns
Returns the computed energies based on the input parameters.

# Notes
Ensure that the dimensions of `weights` and `localEnergies` are compatible.
"""
function getEnergies(weights,localEnergies,PMax;
    Gnp = precomputeNormalizedAccWeight(weights,PMax),)
    
    N = lastindex(localEnergies)
    num = zeros(PMax)
    denom = zeros(PMax)
    for p in 1:PMax
        for n in p+1:N
            num[p] += Gnp[n,p]*localEnergies[n]
            denom[p] += Gnp[n,p]
        end
    end
    return num ./denom
end

function getEnergies(Obs::BasicObserver,PMax;n_equilibration=1)
    @views getEnergies(Obs.TotalWeights[n_equilibration:end],Obs.energies[n_equilibration:end],PMax)
end

function getEnergies(Obs::CombinedObserver,PMax;kwargs...);
    for O in Obs.Observers
        O isa BasicObserver && return getEnergies(O,PMax;kwargs...)
    end
    error("No BasicObserver found in CombinedObserver. Consider explicitly running `getEnergies(weights,localEnergies,PMax)` instead")
end
