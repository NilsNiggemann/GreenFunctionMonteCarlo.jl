function RandomReset!(prob;ratio = 0.7)
    confs = getConfigs(prob)
    for i in eachindex(confs)
        if rand() < ratio
            randomRydbergState!(confs[i], Jij_Mat)
        end
    end
    return prob
end
function getAllConfs(P::ProblemEnsemble)
    confs = Set{eltype(P.problems[begin].Walkers.Configs)}()
    for i in eachindex(P.problems)
        Pconfs = getConfigs(P.problems[i]).x
        
        for j in eachindex(Pconfs)
            for c in Pconfs[j]
                push!(confs, c)
            end
        end
    end
    return confs
end

struct LocalEnergy{HType <: LocalOperator,GWFType <: AbstractGuidingFunction, HT <: AbstractHilbertSpace}
    H::HType
    logpsi::GWFType
    Hilbert::HT
end
function (le::LocalEnergy)(x::AbstractVector)
    return getLocalEnergy(x,le.H, NaiveFunction(le.logpsi), le.Hilbert)
end

function P_exchange(P1::GFMCProblem,x1,P2::GFMCProblem,x2)
    el¹ = LocalEnergy(P1.H, P1.logψ, P1.Hilbert)
    el² = LocalEnergy(P2.H, P2.logψ, P2.Hilbert)
    return P_exchange(el¹(x1), el²(x2), el²(x2), el¹(x1))
end

function P_exchange(e1::LocalEnergy,e2::LocalEnergy,x1,x2)
    delta = (e2(x1) - e1(x1)) + (e1(x2) - e2(x2))
    ratio = delta/abs(e2(x2) + e1(x1))
    return min(1.0, exp(-ratio))
end

function inplaceSwap!(A1,A2)
    for i in eachindex(A1,A2)
        A1[i], A2[i] = A2[i], A1[i]
    end
    return A1, A2
end

function replicaExchange!(P::ProblemEnsemble)
    for i_prob in eachindex(P.problems)[1:end-1]
        P_left = P.problems[i_prob]
        P_right = P.problems[i_prob+1]
        confs_left = getConfigs(P_left).x
        confs_right = getConfigs(P_right).x
        for x1 in confs_left
            for x2 in confs_right
                prob = P_exchange(P_left, x1, P_right, x2)
                if rand() < prob
                    inplaceSwap!(x1, x2)
                end
            end
        end
        prob = P_exchange(P_left, confs_left, P_right, confs_right)
    end
end

function chunk_range(NSteps, n_chunks)
    chunk_size = div(NSteps, n_chunks)
    return [(i-1)*chunk_size+1:min(i*chunk_size, NSteps) for i in 1:n_chunks]
end