using DataStructures
using OrderedCollections
using DataFrames
using HTTP
import Dates
using Printf
using JSON3
using SimpleLog


function dates_to_uts(date)
    floor(Int64, Dates.datetime2unix(date))
end


###> start download
"""
    download_ise_power_data(hdf5_dir, year)

    hdf5_dir : directory where to store downloaded json file
    year : power data as json is downloaded for year and saved to power_<year>.json

    /public_power Public Power
    /cbpf Cross Border Physical Flows
"""
function download_ise_power_data(json_dir, year_, get_)
    query = [("country" => "de"), ("start", @sprintf("%s-01-01T00:00Z",year_)), ("end", @sprintf("%s-12-31T23:59Z", year_))]
    url   = "https://api.energy-charts.info/"*get_
    
    @infoe @sprintf("Download %s: %d", get_, year_)
    resp  = HTTP.get(url, ["Accept" => "application/json"], query = query)

    b     = String(resp.body)
    bb    = JSON3.read(b)
    path  = joinpath(json_dir,@sprintf("%s_%s.json", get_, year_))
    open(path, "w") do out
        JSON3.pretty(out, bb)
    end
    @infoe @sprintf("saved to %s", path)
end

function load_ise_energy_chart_data(json_dir, start_year, end_year, get)
    for year in start_year:end_year
        download_ise_power_data(json_dir, year, get)
    end
end

"""
    download_ise_istalled_power_data()

    downlaod the intsalled data and save in file installed_power.json
"""
function download_ise_istalled_power_data(json_dir)
    query = Dict("country" => "de", "time_step" => "monthly", "installation_decomission" => false)
    url   = "https://api.energy-charts.info/installed_power"
    resp  = HTTP.get(url, ["Accept" => "application/json"], query = query)
    b = String(resp.body)
    bb = JSON3.read(b)
    
    function make_datetime(x)
        my = split(x, ".")
        Dates.DateTime(parse(Int64, my[2]), parse(Int64, my[1]))
    end

    uts = [dates_to_uts(make_datetime(x)) for x in bb["time"]]

    bbb = Dict("unix_seconds" => uts, "production_types" => bb["production_types"])
    for (n,d) in bbb["production_types"]
        @infoe n
    end

    open(joinpath(json_dir, "installed_power.json"), "w") do out
        JSON3.pretty(out, bbb)
    end
end

"""
    Download ise data and store as json in json_dir
"""
function download_ise_data(json_dir, start_year, end_year)
    # cbpf Cross Border Physical Flows
    gets = ["public_power", "total_power", "installed_power", "cbpf"]
    for get in gets
        load_ise_energy_chart_data(json_dir, start_year, end_year, get)
    end
end
###< end download

jsonroot = joinpath(dataroot)
json_dir, start_year, end_year = dataroot, 2025, 2025

download_ise_data(json_dir, start_year, end_year)