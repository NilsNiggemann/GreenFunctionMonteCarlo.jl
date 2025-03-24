function precomputeNormalizedAccWeight(weights,PMax;normalize=true)
    bn = weights
    meanweight = normalize ? Statistics.mean(bn) : 1
    
    Gnp = zeros(length(bn),PMax)

    for n in axes(Gnp,1)
        Gnp[n,1] = 1 #zero projection order
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

function getEnergies(weights,localEnergies,PMax;
    Gnp = precomputeNormalizedAccWeight(weights,PMax)
    )
    
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

function getEnergies(Obs::BasicObserver,PMax);
    return getEnergies(Obs.TotalWeights,Obs.energies,PMax)
end

function getEnergies(Obs::CombinedObserver,PMax);
    for O in Obs.Observers
        O isa BasicObserver && return getEnergies(O,PMax)
    end
    error("No BasicObserver found in CombinedObserver. Consider explicitly running `getEnergies(weights,localEnergies,PMax)` instead")
end
