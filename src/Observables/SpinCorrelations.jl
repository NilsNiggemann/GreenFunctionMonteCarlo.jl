struct SpinCorrelations{T<:Real} <: AbstractObservable
    obsBuffer::Vector{T}
    i_inds::Vector{Int}
    j_inds::Vector{Int}
end
function SpinCorrelations(Nsites::Integer,t=Int8) 
    ij_inds = upTriIter(LinearAlgebra.I(Nsites))
    
    SpinCorrelations(zeros(t,(Nsites * (Nsites + 1))÷2), [i for (i,j) in ij_inds], [j for (i,j) in ij_inds])
end
Base.copy(O::SpinCorrelations) = SpinCorrelations(copy(O.obsBuffer), copy(O.i_inds), copy(O.j_inds))
obs(O::SpinCorrelations) = O.obsBuffer

upTriIter(A::AbstractMatrix) = ((i,j) for i in axes(A,1) for j in i:size(A,2))

function (O::SpinCorrelations)(out,config::AbstractVector)
    ONE = one(eltype(out))
    TWO = 2 * ONE
    sz(n) = ONE - TWO * n

    GreenFunctionMonteCarlo.LoopVectorization.@turbo for index in eachindex(O.i_inds, O.j_inds)
        i = O.i_inds[index]
        j = O.j_inds[index]
        xi = config[i]
        xj = config[j]
        si = sz(xi)
        sj = sz(xj)
        out[index] = si * sj
    end
    return out
end

(O::SpinCorrelations)(out,config::BosonConfig) = O(out, parent(config))

function average_obs_walkers(O::SpinCorrelations,obs_idx::Integer,walker_confs::AbstractMatrix, WalkerPopulations::AbstractVector{<:Integer})
    Obs_num_i_m_b = zero(eltype(obs(O)))

    sz(n) = 1 - 2 * n
    i,j = O.i_inds[obs_idx], O.j_inds[obs_idx]

    LoopVectorization.@turbo for α in axes(walker_confs,1)
        mult = WalkerPopulations[α]
        x_i = walker_confs[α,i]
        x_j = walker_confs[α,j]
        si = sz(x_i)
        sj = sz(x_j)

        Obs_num_i_m_b += si*sj *mult
    end
    return Obs_num_i_m_b
end

""" Helper function to convert the upper triangular buffer of spin correlations into a matrix form for easier analysis. 
    get_matrix_from_tri(els::AbstractVector) -> AbstractMatrix
"""
function get_matrix_from_tri(els::AbstractVector)
    Nsites = Int((sqrt(8length(els)+1)-1)÷2)
    mat = zeros(eltype(els), Nsites, Nsites)
    for (index,(i,j)) in enumerate(upTriIter(mat))
        mat[i,j] = els[index]
        mat[j,i] = els[index]

    end
    return mat
end