using DataStructures
using OrderedCollections
using DataFrames
import Dates
using Printf
using JSON3
using 

function uts_to_dates(uts)
    Dates.unix2datetime(uts)
end

function replace_missing(V::Vector{T}, a::T) where T
    @. ifelse(V == missing, a, V)
end

function extract_start_end_date(hdf5_dir, year)
    path = joinpath(hdf5_dir, @sprintf("public_power_%s.json", year))
    bb = open(path, "r") do io
        JSON3.read(read(io, String))
    end
    uts_key = "unix_seconds"

    date_time = @. unix2datetime(bb[uts_key])
    @info (date_time[1], date_time[end])
    @info get_nb_days(date_time[1], date_time[end])
end

function name_years(name,  start_year, end_year; ext = nothing)
    if ext === nothing
        return @sprintf("%s_%d-%d", name, start_year, end_year)
    end
    @sprintf("%s_%d-%d.%s", name, start_year, end_year, ext)
end

###> 
"""
["Biomass (GW)", "Wind offshore (GW)", "Battery Storage (Capacity) (GWh)", "Wind onshore (GW)", "Battery Storage (Power) (GW)", "months", "Solar (GW)"]
"""
function installed_power_to_hdf5(json_dir, hdf5_dir)
    path = joinpath(json_dir, "installed_power.json")
    IP = open(path, "r") do io
        JSON3.read(read(io, String))
    end

    uts = IP["unix_seconds"]
    @info "installed power:", uts_to_dates(uts[1]), uts_to_dates(uts[end])
    prodtypes = IP["production_types"]
    for (n,d) in prodtypes
        @info n
    end

    ip_keys = ["Wind onshore", "Wind offshore", "Solar",  "Biomass", "Battery Storage (Capacity)", "Battery Storage (Power)"]    
    datasets = Dict()
    for (i,k) in enumerate(ip_keys)
        flag = false
        for p in prodtypes
            if k == p["name"]
                data  = @. ifelse(p["data"]  === nothing, 0.0, p["data"])
                datasets[k] = data
                flag = true
            end
        end
        if !flag
            @info @sprintf("key %s not in data", k)
        end
    end
    datasets["uts"] = Array{Float64}(uts)
    groups = Dict("installed_power" => datasets)
    hdf5_path = joinpath(hdf5_dir, "installed_power.hdf5")

    @info "Save installed_power to", hdf5_path
    save_groups_as_hdf5(hdf5_path, groups, script_dir=false, permute_dims_p = true)
end

function load_installed_power_from_hdf5(hdf5_path)
    groups = load_groups_as_hdf5(hdf5_path, script_dir=false, permute_dims_p = true)
    uts = groups["installed_power"]["uts"]
    data = Dict()
    for (k,d) in groups["installed_power"]
        if k != "uts"
            data[k] = d
        end
    end
    uts, data
end
    
"""
    keys in energy_charts json files

    unix_seconds                     [Unix utstamp]
    Hydro pumped storage consumption [MW]
    Import Balance                   [MW]
    Nuclear                          [MW]
    Hydro Run-of-River               [MW]
    Biomass                          [MW]
    Fossil brown coal / lignite      [MW]
    Fossil coal-derived gas          [MW]
    Fossil hard coal                 [MW]
    Fossil oil                       [MW]
    Fossil gas                       [MW]
    Geothermal                       [MW]
    Hydro water reservoir            [MW]
    Hydro pumped storage             [MW]
    Others                           [MW]
    Waste                            [MW]
    Wind offshore                    [MW]
    Wind onshore                     [MW]
    Solar                            [MW]
    Load                             [MW]
    Residual load                    [MW]
    Renewable share of generation    [%]
    Renewable share of load          [%]

    ise_json_to_dict(uts_key, ise_keys, data_root, hdf5_path1, hdf5_path2)

    read all power data json files, combine them and store in hdf5 format

    uts_key    : key of the time instances column given as Unix utstamp, "xAxisValues (Unix timestamp)"
    ise_keys   : keys for data columns to include
    data_root  : directory where jsons files are
    hdf5_path1 : data matrices per year are stored in this hdf5 file
    hdf5_path2 : concatenated data matrices are stored as singel matrix in this hdf5 file
"""
function ise_json_to_dict(fname, tname, pname, json_dir; start_year::Int64, end_year::Int64)
    datafiles = [@sprintf("%s_%d.json", fname, y) for y in start_year:end_year]

    # dataset_name, data
    datasets = OrderedDict{String, Dict{String, Vector{Float64}}}()
    time_stamps = OrderedDict{String, Vector{Float64}}()

    df = datafiles[1]
    for df in datafiles
        # read and parse json
        path = joinpath(json_dir, df)
        d = open(path, "r") do io
            JSON3.read(io)
        end

        D = OrderedDict{String, Vector{Float64}}()    
        production_types = d[pname]
        for prodtype in production_types
            data = prodtype["data"]
            V = try
                Vector{Float64}(data)
            catch
                V = Vector{Float64}(undef,0)
                #@info ise_keys[i], typeof(d[ise_keys[i]])
                for e in data
                    if typeof(e) == Float64
                        push!(V, e)
                    else
                        push!(V, 0.0)
                    end
                end
                V
            end
            D[prodtype["name"]] = V
        end
        
        fname = splitext(df)[1]
        @info fname
        tt = if eltype(d[tname]) == String
            [parse(Float64, x) for x in d[tname]]
        else
            d[tname]
        end

        time_stamps[fname] = tt
        datasets[fname] = D
    end
    time_stamps, datasets
end

function plot_installed_power(json_dir)
    hp = name_years("ise_installed_power", start_year, end_year, ext = "json")
    d = JSON.parsefile(joinpath(json_dir, hp), dicttype=DataStructures.OrderedDict, inttype=Int64, use_mmap=true)
    ms = d["months"]
    dates = map(x -> Date(parse(Int, x[2]), parse(Int,x[1]), 15), split.(ms, "."))
    Date(2012,2,29)

    bm  = d["Biomass (GW)"]                     #
    wof = d["Wind offshore (GW)"]               #
    Won = d["Wind onshore (GW)"]                #
    PV  = d["Solar (GW)"]                       #
    BP  = d["Battery Storage (Power) (GW)"]     #
    BC  = d["Battery Storage (Capacity) (GWh)"] #

    bm  = @. ifelse(bm  === nothing, 0.0, bm )
    wof = @. ifelse(wof === nothing, 0.0, wof)
    Won = @. ifelse(Won === nothing, 0.0, Won)
    PV  = @. ifelse(PV  === nothing, 0.0, PV )
    BP  = @. ifelse(BP  === nothing, 0.0, BP )
    BC  = @. ifelse(BC  === nothing, 0.0, BC )
    plt.plot(dates, bm , label="Bio")
    plt.plot(dates, wof, label="Wof")
    plt.plot(dates, Won, label="Won")
    plt.plot(dates, PV , label="PV")
    plt.plot(dates, BP , label="BP")
    plt.plot(dates, BC , label="BC")
    plt.legend()
end

function save_ise_data_to_hdf5(hdf5_dir, hdf5_name, time_stamps, datasets, start_year, end_year)
    groups = Dict{String, Dict{String, Vector{Float64}}}()
    lk = []
    ld = []
    prodtypes = OrderedSet{String}()
    for (fname, D) in datasets
        group = Dict{String, Vector{Float64}}()
        group["time_stamp"] = time_stamps[fname]
        push!(lk, length(keys(D)))
        for (prodtype, data) in D
            push!(prodtypes, prodtype)
            group[prodtype] = data
            push!(ld, length(data))
        end
        groups[fname] = group
    end
    hdf5_path1 = joinpath(hdf5_dir, @sprintf("%s_per_year.hdf5", hdf5_name))
    save_groups_as_hdf5(hdf5_path1, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)

    prodtypes = collect(prodtypes)
    dpr = OrderedDict()
    lpr = []
    for (i, pr) in enumerate(prodtypes)
        dpr[pr] = i
        push!(lpr, (i, pr))
    end
    sort!(lpr, by=x -> x[1])
    prs = [a[2] for a in lpr]
    
    matmat = []
    for (fname, D) in datasets
        ld = length(time_stamps[fname])
        mat = zeros(Float64, ld, maximum(lk)+1)
        mat[:,1] = time_stamps[fname]
        for (prodtype, data) in D
            i = dpr[prodtype]
            #@info fname, i, prodtype
            mat[:,i+1] = data
        end
        #@info size(mat)
        push!(matmat, mat)
    end
    DD = reduce(vcat, matmat)
    groups = Dict(@sprintf("%s_2016_2024",hdf5_name) => Dict("data" => DD, "colnames" => prs))

    hp =  @sprintf("%s_%d-%d.%s", hdf5_name, start_year, end_year, ext="hdf5")
    hdf5_path2 = joinpath(hdf5_dir, hp)
    save_groups_as_hdf5(hdf5_path2, groups; permute_dims_p=true, extension=".hdf5", script_dir=false)

    hdf5_path2
end

function load_ise_data_from_hdf5(hdf5_path)
    groups = load_groups_as_hdf5(hdf5_path; permute_dims_p=true)
    grname = collect(keys(groups))[1]
    data = groups[grname]["data"]
    prodtypes = groups[grname]["colnames"]
    data, prodtypes
end

"""
    downlaod data from enerycharts and convert to hdf5
"""
function ise_json_to_hdf5(json_dir, hdf5_dir, start_year, end_year)
    time_stamps, datasets = ise_json_to_dict("public_power", "unix_seconds", "production_types", json_dir, start_year=start_year, end_year=end_year);
    public_power_hdf5_path = save_ise_data_to_hdf5(hdf5_dir, "public_power", time_stamps, datasets, start_year, end_year)
    
    time_stamps, datasets = ise_json_to_dict("total_power", "unix_seconds","production_types", json_dir, start_year=start_year, end_year=end_year);
    total_power_hdf5_path = save_ise_data_to_hdf5(hdf5_dir, "total_power", time_stamps, datasets, start_year, end_year)

    time_stamps, datasets = ise_json_to_dict("cbpf", "unix_seconds", "countries", json_dir, start_year=start_year, end_year=end_year);
    cbpf_hdf5_path = save_ise_data_to_hdf5(hdf5_dir,"cbpf", time_stamps, datasets, start_year, end_year)

    installed_power_to_hdf5(json_dir, hdf5_dir)
end

function get_paths()
    # project root directory
    sisremo_dir = dirname(dirname(@__DIR__))
    hdf5_dir = joinpath(sisremo_dir, "data")
    json_dir = joinpath(hdf5_dir, "json_downloads")

    mkpath(hdf5_dir)
    mkpath(json_dir)
    json_dir, hdf5_dir
end

function json_to_hdf5()
    start_year = 2016
    end_year = 2025
    json_dir, hdf5_dir = get_paths()

    #download_ise_data(json_dir, start_year, end_year)
    #download_ise_istalled_power_data(json_dir)

    ise_json_to_hdf5(json_dir, hdf5_dir, start_year, end_year)
end
