
struct ArrayWalkerEnsemble{T,N,BuffType} <: AbstractWalkerEnsemble 
    Configs::AbstractArray{T,N}
    weights::AbstractVector{T}
    buffers::AbstractVector{BuffType}
end

parent(X::ArrayWalkerEnsemble) = X.Configs

getConfigs(X::ArrayWalkerEnsemble{T,N}) where {T,N} = eachslice(X,dims=N)