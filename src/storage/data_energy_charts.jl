using DataStructures
using HTTP
using Dates
using Printf
using Common # json io, @ionfoe

import PyPlot as plt
plt.pygui(true)
plt.pygui(:qt5)

# project root directory
sisremo_dir = dirname(dirname(@__DIR__))
get_data_dir() = joinpath(sisremo_dir, "data")
mkpath(get_data_dir())

"""
    download_ise_power_data(data_dir, year)

    data_dir : directory where to store downloaded json file
    year : power data as json is downloaded for year and saved to power_<year>.json
"""
function download_ise_power_data(data_dir, year)
    query = [("country" => "DE"), ("start", @sprintf("%s-01-01T00:00Z",year)), ("end", @sprintf("%s-12-31T23:59Z", year))]
    url   = "https://api.energy-charts.info/power"
    resp  = HTTP.get(url, ["Accept" => "application/json"], query = query)
    str   = String(resp.body)
    path  = joinpath(data_dir, @sprintf("power_%s.json", year))
    # write String str to path
    open(path, "w") do out
        write(out, str)
    end
end
function load_ise_energy_chart_data(data_dir, start_year, end_year)
    for year in start_year:end_year
        download_ise_power_data(data_dir, year)
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
#extract_start_end_date(get_data_dir(), year)

"""
    download_ise_istalled_power_data()

    downlaod the intsalled data and save in file installed_power.json
"""
function download_ise_istalled_power_data(data_dir)
    query = Dict("country" => "de", "time_step" => "monthly", "installation_decomission" => false)
    url   = "https://api.energy-charts.info/installed_power"
    resp  = HTTP.get(url, ["Accept" => "application/json"], query = query)
    b     = String(resp.body)
    open(joinpath(data_dir, "installed_power.json"), "w") do out
        write(out, b)
    end
end
#download_ise_istalled_power_data(get_data_dir())

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

    dates  = @. year_month_to_date(IP["months"])
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
    println(dates[1], " ", dates[end])
    months = @. dates_to_uts(dates)

    ip_keys = ["Wind onshore (GW)", "Wind offshore (GW)", "Solar (GW)",  "Biomass (GW)", "Battery Storage (Capacity) (GWh)", "Battery Storage (Power) (GW)"]
    M = Matrix{Float64}(undef, length(months), length(ip_keys)+1)

    M[:,1] = months
    for (i,k) in enumerate(ip_keys)
        data  = @. ifelse(IP[k]  === nothing, 0.0, IP[k])
        M[:,i+1] = data[i1:i2]
    end
    hp = name_years("ise_installed_power", start_year, end_year, ext = "hdf5")
    hdf5_path = joinpath(data_dir, hp)
    save_array_as_hdf5(hdf5_path, M, group = "ise_installed_power", dataset = "ise_installed_power", script_dir=false, colm_to_rowm_p = true)
end
#installed_power_to_hdf5(get_data_dir(), 2016, 2022)

"""
    keys in energy_charts json files

    xAxisValues                      [Unix utstamp]
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

    ise_json_to_hdf5(uts_key, ise_keys, data_root, hdf5_path1, hdf5_path2)

    read all power data json files, combine them and store in hdf5 format

    uts_key    : key of the time instances column given as Unix utstamp, "xAxisValues (Unix timestamp)"
    ise_keys   : keys for data columns to include
    data_root  : directory where jsons files are
    hdf5_path1 : data matrices per year are stored in this hdf5 file
    hdf5_path2 : concatenated data matrices are stored as singel matrix in this hdf5 file
"""
function ise_json_to_hdf5(uts_key, ise_keys, data_root, datafiles, hdf5_path1, hdf5_path2)
    # number of columns
    nk = length(ise_keys)
    datasets = []
    # destination of data
    D = Vector{Array{Float64,2}}(undef, 0)

    ntot = 0
    for df in datafiles
        @infoe df
        # read and parse json
        path = joinpath(data_root, df)
        d = from_json(path)
        # d ise_keys to include
        @infoe keys(d)

        # length of data vector
        n = length(d[uts_key])
        ntot += n
        # create Matrix
        M = Matrix{Float64}(undef, n, nk+1)
        # store Unix utstamps
        M[:,1] = d[uts_key]
        # get data and store in Matrix M
        for (i,k) in enumerate(ise_keys)
            if haskey(d, ise_keys[i])
                V = d[ise_keys[i]]
                V = replace_missing(V,  toreplace = nothing)
                M[:,i+1] = V
            else
                @infoe "no", ise_keys[i]
                M[:,i+1] .= 0.0
            end
        end
        # push Matrix into Vector
        push!(D, M)
        push!(datasets, split(splitext(df)[1], "_")[2])
    end
    # save al matrices in a single hdf5file
    save_arrays_as_hdf5(hdf5_path1, D, group = "ise_power", datasets = datasets, script_dir=false, colm_to_rowm_p = true)

    # concatenate all matrices
    DD = reduce(vcat, D)

    # save the single data matrix DD
    save_array_as_hdf5(hdf5_path2, DD, group = "ise_power", dataset = "ise_power", script_dir=false, colm_to_rowm_p = true)
end

"""
    run_ise_json_to_hdf5(download_json_p::Bool, start_year::Int64, end_year::Int64)

    download_json_p : if true downlaod energy_chart data start_year:end_year

"""
function run_ise_json_to_hdf5(data_dir, download_json_p::Bool, start_year::Int64, end_year::Int64)
    ise_keys_all = [
        "Load (MW)",                                #  1
        "Residual load (MW)",                       #  2
        "Wind offshore (MW)",                       #  3
        "Wind onshore (MW)",                        #  4
        "Solar (MW)",                               #  5
        "Biomass (MW)",                             #  6
        "Hydro pumped storage consumption (MW)",    #  7
        "Hydro water reservoir (MW)",               #  8
        "Hydro pumped storage (MW)",                #  9
        "Hydro Run-of-River (MW)",                  # 10
        "Import Balance (MW)",                      # 11
        "Geothermal (MW)",                          # 12
        "Waste (MW)",                               # 13
        "Nuclear (MW)",                             # 14
        "Fossil brown coal / lignite (MW)",         # 15
        "Fossil coal-derived gas (MW)",             # 16
        "Fossil hard coal (MW)",                    # 17
        "Fossil oil (MW)",                          # 18
        "Fossil gas (MW)",                          # 19
        "Others (MW)",                              # 20
        "Renewable share of generation (%)",        # 21
        "Renewable share of load (%)"]              # 22

    if download_json_p
        load_ise_energy_chart_data(data_dir, start_year, end_year)
    end

    uts_key    = "xAxisValues (Unix timestamp)"
    datafiles  = [@sprintf("power_%d.json", y) for y in start_year:end_year]

    hdf5_path1 = joinpath(data_dir, "ise_power_all.hdf5")
    hp =  @sprintf("ise_power_all_%d-%d.%s", start_year, end_year, ext="hdf5")
    hdf5_path2 = joinpath(data_dir,hp)

    ise_json_to_hdf5(uts_key, ise_keys_all, data_dir, datafiles, hdf5_path1, hdf5_path2)
end
#run_ise_json_to_hdf5(get_data_dir(), false, 2016, 2022)

"""
    load ise energy charts data stored in hdf5 file
"""
function load_ise_as_hdf5(data_dir, start_year, end_year)
    hp = name_years("ise_power_all", start_year, end_year, ext = "hdf5")
    hdf5_path = joinpath(data_dir, hp)
    load_array_as_hdf5(hdf5_path, group = "ise_power", dataset = "ise_power", script_dir=false, colm_to_rowm_p = true)
end

function load_ise_installed_power(data_dir, start_year, end_year)
    hp = @sprintf("ise_installed_power_%s-%s.hdf5", start_year, end_year)
    hdf5_path = joinpath(data_dir, hp)
    load_array_as_hdf5(hdf5_path, group = "ise_installed_power", dataset = "ise_installed_power", script_dir=false, colm_to_rowm_p = true)
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
