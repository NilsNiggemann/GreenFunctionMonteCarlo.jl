function createMMapArray(file::HDF5.File,datasetname::String,type,dims)
    SaveConfigs_dset = HDF5.create_dataset(file,datasetname,HDF5.datatype(type),HDF5.dataspace(dims);alloc_time = HDF5.H5D_ALLOC_TIME_EARLY)
    @assert HDF5.ismmappable(SaveConfigs_dset) "Dataset is not mappable for given type $(eltype(SaveConfigs_dset))"
    return HDF5.readmmap(SaveConfigs_dset)
end
function maybe_MMap_array(filename::AbstractString,datasetname::String,type,dims)
    HDF5.h5open(filename,"cw") do file
        return createMMapArray(file,datasetname,type,dims)
    end
end

function maybe_MMap_array(filename::Nothing,datasetname::String,type,dims)
    return zeros(type,dims)
end

function readMMapArray(filename::AbstractString,datasetname::String)
    HDF5.h5open(filename,"r") do file
        SaveConfigs_dset = file[datasetname]
        return HDF5.readmmap(SaveConfigs_dset)
    end
end

strd(x,args...;kwargs...) = string(round(x,args...;digits = 3,kwargs...))

"""
    sample_fast(rng, weights) -> Int

Allocation-free categorical sampling. Returns an index sampled proportionally
to `weights`. Functionally equivalent to
`StatsBase.sample(rng, StatsBase.Weights(weights))` but avoids heap-allocating
a `Weights` wrapper object on every call.

This matters especially inside the innermost Monte Carlo loop
(`performMarkovStep!`) where the per-call allocation would otherwise
accumulate into significant GC pressure and hurt multi-threaded scaling.
"""
@inline function sample_fast(rng::Random.AbstractRNG, weights::AbstractVector)
    total = sum(weights)
    u = rand(rng) * total
    cumulative = zero(total)
    @inbounds for i in eachindex(weights)
        cumulative += weights[i]
        u <= cumulative && return i
    end
    return lastindex(weights)  # fallback for floating-point rounding
end


