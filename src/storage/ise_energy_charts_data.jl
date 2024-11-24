using DataStructures
using OrderedCollections
using DataFrames
using HTTP
import Dates
using Printf
using JSON3

using BaseUtils

# project root directory
sisremo_dir = dirname(dirname(@__DIR__))
get_data_dir() = joinpath(sisremo_dir, "data")
data_dir = get_data_dir()
mkpath(get_data_dir())

"""
    download_ise_power_data(data_dir, year)

    data_dir : directory where to store downloaded json file
    year : power data as json is downloaded for year and saved to power_<year>.json
"""
function download_ise_power_data(data_dir, year, get)
    query = [("country" => "de"), ("start", @sprintf("%s-01-01T00:00Z",year)), ("end", @sprintf("%s-12-31T23:59Z", year))]
    url   = "https://api.energy-charts.info/"*get
    @infoe @sprintf("Download %s: %d", get, year)
    resp  = HTTP.get(url, ["Accept" => "application/json"], query = query)
    str   = String(resp.body)
    path  = joinpath(data_dir, @sprintf("%s_%s.json", get, year))
    # write String str to path
    open(path, "w") do out
        write(out, str)
    end
    @infoe @sprintf("saved to %s", path)
end
function load_ise_energy_chart_data(data_dir, start_year, end_year, get)
    for year in start_year:end_year
        download_ise_power_data(data_dir, year, get)
    end
end

function extract_start_end_date(data_dir, year)
    uts_key = "xAxisValues (Unix timestamp)"
    path = joinpath(data_dir, @sprintf("power_%s.json", year))

    b = load_json(path)

    date_time = @. unix2datetime(b[uts_key])
    @info (date_time[1], date_time[end])
    @info get_nb_days(date_time[1], date_time[end])
end

"""
    download_ise_istalled_power_data()

    downlaod the intsalled data and save in file installed_power.json
"""
function download_ise_istalled_power_data(data_dir)
    query = Dict("country" => "de", "time_step" => "monthly", "installation_decomission" => false)
    url   = "https://api.energy-charts.info/installed_power"
    resp  = HTTP.get(url, ["Accept" => "application/json"], query = query)
    b     = String(resp.body)
    bb = JSON3.read(b)
    for (n,d) in bb["production_types"]
        @infoe n
    end
    open(joinpath(data_dir, "installed_power.json"), "w") do out
        #write(out, b)
        JSON3.pretty(out, b)
    end
end

function year_month_to_date(ym)
    m,y = split(ym, ".")
    date = DateTime(parse(Int64, y), parse(Int64, m))
end
function dates_to_uts(date)
    floor(Int64, datetime2unix(date))
end

function name_years(name,  start_year, end_year; ext = nothing)
    if ext === nothing
        return @sprintf("%s_%d-%d", name, start_year, end_year)
    end
    @sprintf("%s_%d-%d.%s", name, start_year, end_year, ext)
end

"""
["Biomass (GW)", "Wind offshore (GW)", "Battery Storage (Capacity) (GWh)", "Wind onshore (GW)", "Battery Storage (Power) (GW)", "months", "Solar (GW)"]
"""
function installed_power_to_hdf5(data_dir, start_year, end_year)
    path = joinpath(data_dir, "installed_power.json")
    IP = from_json(path)
    for (n,d) in IP["production_types"]
        @infoe n
    end

    dates  = @. year_month_to_date(IP["time"])
    i1, i2 = 0, 0
    for (i,d) in enumerate(dates)
        if d < DateTime(start_year-1,12,31)
            i1 = i
        end
        if d < DateTime(end_year+1,1,1)
            i2 = i
        end
    end
    i1 += 1
    dates = dates[i1:i2]
    #println(dates[1], " ", dates[end])
    months = @. dates_to_uts(dates)

    ip_keys = ["Wind onshore", "Wind offshore", "Solar",  "Biomass", "Battery Storage (Capacity)", "Battery Storage (Power)"]
    #ip_keys = ["Wind onshore (GW)", "Wind offshore (GW)", "Solar (GW)",  "Biomass (GW)", "Battery Storage (Capacity) (GWh)", "Battery Storage (Power) (GW)"]
    M = Matrix{Float64}(undef, length(months), length(ip_keys)+1)

    # [Dict"("name"=>, "data"=>)]
    PT = IP["production_types"]

    M[:,1] = months
    for (i,k) in enumerate(ip_keys)
        flag = false
        for p in PT
            if k == p["name"]
                data  = @. ifelse(p["data"]  === nothing, 0.0, p["data"])
                M[:,i+1] = data[i1:i2]
                flag = true
            end
        end
        if !flag
            @infoe @sprintf("key %s not in data", k)
        end
    end
    hp = name_years("ise_installed_power", start_year, end_year, ext = "hdf5")
    hdf5_path = joinpath(data_dir, hp)
    save_array_as_hdf5(hdf5_path, M, group_name = "ise_installed_power", dataset_name = "ise_installed_power", script_dir=false, colm_to_rowm_p = true)
end

function replace_missing(V::Vector{T}, a::T) where T
    @. ifelse(V == missing, a, V)
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
function ise_json_to_dict(fname, tname, pname, data_dir; start_year::Int64, end_year::Int64)
    datafiles  = [@sprintf("%s_%d.json", fname, y) for y in start_year:end_year]

    # dataset_name, data
    datasets = OrderedDict{String, Dict{String, Vector{Float64}}}()
    time_stamps = OrderedDict{String, Vector{Float64}}()

    df = datafiles[1]
    for df in datafiles
        # read and parse json
        path = joinpath(data_dir, df)
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
                #@infoe ise_keys[i], typeof(d[ise_keys[i]])
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
        @infoe fname
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


"""
    load ise energy charts data stored in hdf5 file
"""
function load_ise_as_hdf5(data_dir, start_year, end_year)
    hp = name_years("ise_power_all", start_year, end_year, ext = "hdf5")
    hdf5_path = joinpath(data_dir, hp)
    load_array_as_hdf5(hdf5_path, group_name = "ise_power", dataset_name = "ise_power", script_dir=false, permute_dims_p = true)
end

function load_ise_installed_power(data_dir, start_year, end_year)
    hp = @sprintf("ise_installed_power_%s-%s.hdf5", start_year, end_year)
    hdf5_path = joinpath(data_dir, hp)
    load_array_as_hdf5(hdf5_path, group_name = "ise_installed_power", dataset_name = "ise_installed_power", script_dir=false, permute_dims_p = true)
end

function plot_power()
    uts = floor.(Int, D[:,1])
    dates = unix2datetime.(floor.(Int, D[:,1]))

    averaging_hours = 7*24

    EE = D[:,4] + D[:,5] + D[:,6]# + D[:,7]
    EEav, datesav = averaging_mean(EE, dates, averaging_hours)
    Lav, datesav = averaging_mean(D[:,2], dates, averaging_hours)

    plt.plot(datesav, Lav, label="Load")
    #plt.plot(D[:,1], D[:,3], label="LoadRes")
    plt.plot(datesav, EEav, label="EE")
    #plt.plot(D[:,1], D[:,4], label="Woff")
    #plt.plot(D[:,1], D[:,5], label="Won")
    plt.legend()
end
#plot_power()

function installed_power(data_dir)
    jp = name_years("ise_installed_power", start_year, end_year, ext = "json")
    d = JSON.parsefile(joinpath(data_dir, hp), dicttype=DataStructures.OrderedDict, inttype=Int64, use_mmap=true)
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

function save_ise_data_to_hdf5(hdf5_name, time_stamps, datasets)
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
    hdf5_path1 = joinpath(data_dir, @sprintf("%s_per_year.hdf5", hdf5_name))
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
            #@infoe fname, i, prodtype
            mat[:,i+1] = data
        end
        @infoe size(mat)
        push!(matmat, mat)
    end
    DD = reduce(vcat, matmat)
    groups = Dict(@sprintf("%s_2016_2024",hdf5_name) => Dict("data" => DD, "colnames" => prs))

    hp =  @sprintf("%s_%d-%d.%s", hdf5_name, start_year, end_year, ext="hdf5")
    hdf5_path2 = joinpath(data_dir,hp)
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
function download_ise_data(;download_data = false, start_year=2016, end_year = 2024)
    if download_data
        gets = ["public_power", "total_power", "installed_power", "cbpf"]
        for get in gets
            load_ise_energy_chart_data(data_dir, start_year, end_year, get)
        end
    end

    start_year, end_year = 2016, 2024
    time_stamps, datasets = ise_json_to_dict("public_power", "unix_seconds", "production_types", get_data_dir(), start_year=2016, end_year=2024);
    public_power_hdf5_path = save_ise_data_to_hdf5("public_power", time_stamps, datasets)
    
    time_stamps, datasets = ise_json_to_dict("total_power", "unix_seconds","production_types",get_data_dir(), start_year=2016, end_year=2024);
    total_power_hdf5_path = save_ise_data_to_hdf5("total_power", time_stamps, datasets)

    time_stamps, datasets = ise_json_to_dict("installed_power", "time", "production_types", get_data_dir(), start_year=2016, end_year=2024);
    installed_power_hdf5_path = save_ise_data_to_hdf5("installed_power", time_stamps, datasets)

    time_stamps, datasets = ise_json_to_dict("cbpf", "unix_seconds", "countries", get_data_dir(), start_year=2016, end_year=2024);
    cbpf_hdf5_path = save_ise_data_to_hdf5("cbpf", time_stamps, datasets)

    data, prodtypes = load_ise_data_from_hdf5(hdf5_path)

end
download_ise_data(download_data = false, start_year = 2016, end_year =2024)

