using HDF5

function colm_to_rowm(A)
    permutedims(A, reverse(1:ndims(A)))
end

function save_array_as_hdf5(hdf5_path::String, C::Array{T,N}; group = "A", dataset = "A", script_dir=false, colm_to_rowm_p = false) where {T, N}
    if script_dir
        hdf5_path = joinpath(@__DIR__, hdf5_path)
    end
    local gid
    h5open(hdf5_path, "w") do fid
        try
            gid = g_create(fid, group)
        catch
            gid = create_group(fid, group)
        end
        if colm_to_rowm_p
            gid[dataset] = colm_to_rowm(C)
        else
            gid[dataset] = C
        end
    end
end

function load_array_as_hdf5(hdf5_path; group = "A", dataset = "A", script_dir=false, colm_to_rowm_p = false)
    if script_dir
        hdf5_path = joinpath(@__DIR__, hdf5_path)
    end
    local A
    h5open(hdf5_path, "r") do fid
        gid = fid[group]
        A = read(gid[dataset])
    end
    if colm_to_rowm_p
        try 
            return colm_to_rowm(A)
        catch
            return A
        end
    end
    A
end

function save_arrays_as_hdf5(hdf5_path::String, C::Vector{Array{T,N}}; group = "A", datasets = Vector{String}, script_dir=false, colm_to_rowm_p = false) where {T, N}
    if script_dir
        hdf5_path = joinpath(@__DIR__, hdf5_path)
    end
    local gid
    h5open(hdf5_path, "w") do fid
        try
            gid = g_create(fid, group)
        catch
            gid = create_group(fid, group)
        end
        for (i, M) in enumerate(C)
            if colm_to_rowm_p
                gid[datasets[i]] = colm_to_rowm(M)
            else
                gid[datasets[i]] = M
            end
        end
    end
end
