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

function log_replica_exchange_weight(P1::GFMCProblem,P2::GFMCProblem,x1,x2)
    el¹ = LocalEnergy(P1.H, P1.logψ, P1.Hilbert)
    el² = LocalEnergy(P2.H, P2.logψ, P2.Hilbert)
    return log_replica_exchange_weight(el¹(x1), el²(x2), el²(x2), el¹(x1))
end

function log_replica_exchange_weight(e1::LocalEnergy,e2::LocalEnergy,x1,x2)
    delta = (e2(x1) - e1(x1)) + (e1(x2) - e2(x2))
    ratio = -delta/abs(e2(x2) + e1(x1))
    # @info "" delta ratio e1(x1) e2(x1) e1(x2) e2(x2)
    return ratio
end

# function log_replica_exchange_weight(e1_x1::Real,e2_x1::Real,e2_x2::Real,e1_x2::Real)
#     delta = (e2_x1 - e1_x1) + (e1_x2 - e2_x2)
#     ratio = -delta/abs(e2_x2 + e1_x1)
#     return ratio
# end

function inplaceSwap!(A1,A2)
    for i in eachindex(A1,A2)
        A1[i], A2[i] = A2[i], A1[i]
    end
    return A1, A2
end

function get_replica_energies(P::ProblemEnsemble)
    Nproblems = length(P.problems)
    Energies = zeros(Nproblems,Nproblems)
    for i_prob in eachindex(P.problems)
        P_i = P.problems[i_prob]
        E_i = LocalEnergy(P_i.H, P_i.logψ, P_i.Hilbert)

        for j_prob in eachindex(P.problems)
            P_j = P.problems[j_prob]
            confs = getConfigs(P_j)
            for x in confs
                Energies[x,j_prob,i_prob] = E_i(x)
            end
        end
    end
    return Energies
end

function get_swap_weights(P::ProblemEnsemble)
    Nw = length(eachindex(P.problems[1].Walkers))
    Nproblems = length(P.problems)
    @assert all(x->length(eachindex(x.Walkers)) == Nw, P.problems) "All problems must have the same number of walkers"
    swap_weights = zeros(Nw,Nproblems,Nw,Nproblems)

    for i_prob in eachindex(P.problems)
        P_i = P.problems[i_prob]
        for j_prob in eachindex(P.problems)
            P_j = P.problems[j_prob]
            if i_prob == j_prob
                for α in 1:Nw
                    swap_weights[α,i_prob, α ,j_prob] = 1.0
                end
                continue
            end
            
            for α in 1:Nw
                for β in 1:Nw
                    e1 = LocalEnergy(P_i.H, P_i.logψ, P_i.Hilbert)
                    e2 = LocalEnergy(P_j.H, P_j.logψ, P_j.Hilbert)
                    x1 = getConfigs(P_i)[α]
                    x2 = getConfigs(P_j)[β]
                    # e1x1 = Energies[α,j_prob,i_prob]
                    # e1x2 = Energies[β,i_prob,j_prob]
                    # e2x1 = Energies[β,i_prob,j_prob]
                    # e2x2 = Energies[α,j_prob,i_prob]


                    # @assert e1x1 = e1(x1)
                    # @assert e2x1 = e2(x1)
                    # @assert e1x2 = e1(x2)
                    # @assert e2x2 = e2(x2)

                    swap_weights[α,i_prob,β,j_prob] = exp(log_replica_exchange_weight(e1,e2,x1,x2))
                end
            end
        end
    end
    return swap_weights
end

function replicaExchange!(P::ProblemEnsemble;rng = Random.default_rng())
    Nw = length(eachindex(P.problems[1].Walkers))
    @assert all(x->length(eachindex(x.Walkers)) == Nw, P.problems) "All problems must have the same number of walkers"

    Nproblems = length(P.problems)

    swap_weights = get_swap_weights(P)

    for i_prob in eachindex(P.problems)
        P_i = P.problems[i_prob]
        confs = getConfigs(P_i)
        for α in 1:Nw
            Weights = @view swap_weights[α,i_prob,:,:]
            Weights_linear = reshape(Weights, length(Weights))
            new_idx = StatsBase.sample(rng,StatsBase.Weights(Weights_linear))
            
            β,j_prob = Tuple(CartesianIndices(Weights)[new_idx])
            P_j = P.problems[j_prob]
            x1 = confs[α]
            x2 = getConfigs(P_j)[β]

            inplaceSwap!(x1,x2)
        end
    end
    return P
end

function chunk_range(NSteps, n_chunks)
    chunk_size = div(NSteps, n_chunks)
    return [(i-1)*chunk_size+1:min(i*chunk_size, NSteps) for i in 1:n_chunks]
end