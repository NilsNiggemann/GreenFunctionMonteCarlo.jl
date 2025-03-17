function createMMapArray(file::HDF5.File,datasetname::String,type,dims)
    SaveConfigs_dset = create_dataset(file,datasetname,datatype(type),dataspace(dims);alloc_time = HDF5.H5D_ALLOC_TIME_EARLY)
    @assert HDF5.ismmappable(SaveConfigs_dset) "Dataset is not mappable for given type $(eltype(InitConfig))"
    return HDF5.readmmap(SaveConfigs_dset)
end

function readMMapArray(filename::AbstractString,datasetname::String)
    h5open(filename,"r") do file
        SaveConfigs_dset = file[datasetname]
        return HDF5.readmmap(SaveConfigs_dset)
    end
end