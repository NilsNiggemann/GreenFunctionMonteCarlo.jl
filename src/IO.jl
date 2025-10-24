"""
    save_params_dict(file::Union{HDF5.File, HDF5.Group}, name::AbstractString, data)

Save data to an HDF5 file or group.

# Arguments
- `file`: The HDF5 file or group to save to
- `name`: The name of the dataset to create
- `data`: The data to save
"""
function save_params_dict(file::Union{HDF5.File, HDF5.Group}, name, data)
    try
        file[name] = data
    catch e
        @warn "Failed to save $name to HDF5: $e"
    end
end

"""
    save_params_dict(file::Union{HDF5.File, HDF5.Group}, name::AbstractString, data::Dict)

Save a dictionary to an HDF5 file or group. Creates a group with each key-value pair saved recursively.

# Arguments
- `file`: The HDF5 file or group to save to
- `name`: The name of the group to create
- `data`: The dictionary to save
"""
function save_params_dict(file::Union{HDF5.File, HDF5.Group}, name, data::Dict)
    # Create a group for the dictionary
    if isnothing(name) || name == ""
        group = file
    elseif name isa String
        group = HDF5.create_group(file, name)
    end

    # Recursively save each key-value pair
    for (key, value) in data
        save_params_dict(group, string(key), value)
    end
end
save_params_dict(file::Union{HDF5.File, HDF5.Group}, name::AbstractString, data::NamedTuple) = save_params_dict(file, name, Dict(x=>y for (x,y) in pairs(data)))
"""
    save_params_dict(filename::AbstractString, data::Dict; mode::AbstractString="w")

Save a dictionary to an HDF5 file.

# Arguments
- `filename`: The name of the file to save to
- `data`: The dictionary to save
- `mode`: The mode to open the file in ("w" = write/overwrite, "a" = append)
"""
function save_params_dict(filename::AbstractString, data; mode::AbstractString="cw")
    HDF5.h5open(filename, mode) do file
        save_params_dict(file, nothing,data)
    end
end
