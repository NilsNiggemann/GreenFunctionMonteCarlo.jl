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