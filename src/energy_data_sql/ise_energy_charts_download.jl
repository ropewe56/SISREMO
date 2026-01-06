include("include_energy_data.jl")

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
    
    @info @sprintf("Download %s: %d", get_, year_)
    resp  = HTTP.get(url, ["Accept" => "application/json"], query = query)

    b     = String(resp.body)
    bb    = JSON3.read(b)
    path  = joinpath(json_dir,@sprintf("%s_%s.json", get_, year_))
    open(path, "w") do out
        JSON3.pretty(out, bb)
    end
    @info @sprintf("saved to %s", path)
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
        @info n
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

function json_to_dataframes(jsonpaths)

    function normalize_keys(din)
        dout = Dict()
        pk = collect(keys(din))
        for k in pk
            kn = replace(k, " " => "_")
            kn = replace(kn, "/" => "")
            kn = replace(kn, "-" => "_")
            kn = replace(kn, "__" => "_")
            dout[kn] = din[k]
        end
        dout
    end

    alldf = Dict()
    for (table, jsonpath) in jsonpaths

        cont = open(jsonpath, "r") do io
            JSON3.read(io)
        end

        tk = if haskey(cont, :unix_seconds)
            :unix_seconds
        elseif haskey(cont, :time)
            :time
        end
        
        mk = if haskey(cont, :production_types)
            :production_types
        elseif haskey(cont, :countries)
            :countries
        end

        uts = cont[tk]
        din = Dict()
        for pt in cont[mk]
            data = @. ifelse(pt[:data] === nothing, 0.0, pt[:data])
            din[pt[:name]] = Float64.(data)
        end
        dout = normalize_keys(din)

        colnames = Vector{String}(undef,0)
        values   = []
        push!(colnames, string(tk))
        push!(values, uts)
        for (k,v) in dout
            push!(colnames, k)
            push!(values, v)
        end
        df = DataFrame(colnames .=> values)
        
        alldf[table] = df
    end
    alldf
end


dataroot = joinpath(dirname(dirname(@__DIR__)), "data")
jsonroot = joinpath(dataroot, "json")
start_year, end_year = dataroot, 2025, 2025

download_ise_data(jsonroot, start_year, end_year)

jsonpaths= [("cbpf", joinpath(jsonroot, "cbpf_2025.json")),
            ("installed_power", joinpath(jsonroot, "installed_power_2025.json")),
            ("public_power", joinpath(jsonroot, "public_power_2025.json")),
            ("total_power", joinpath(jsonroot, "total_power_2025.json"))
            ]

alldf = json_to_dataframes(jsonpaths)


dbname = joinpath(dataroot, "ise_data.db")
db = SQLite.DB(dbname)

tables = ["cbpf", "public_power", "total_power"]

cmd = "SELECT unix_seconds FROM cbpf;"
df_uts = DBInterface.execute(db, cmd) |> DataFrame
last_uts = maximum(df_uts[!,:unix_seconds])

df_ = filter(row -> row[:unix_seconds] > last_uts, alldf["cbpf"])
alldf["cbpf"][!,:unix_seconds]


table = SQLite.tables(db)[1]
for table in SQLite.tables(db)
    table_name = table.name
    colnames = tables_schema = table.schema.names
    coltypes = tables_schema = table.schema.types
    
    if haskey(alldf, table_name)
        dfnames = Symbol.(names(alldf[table_name]))
        for colname in colnames 
            println(colname in dfnames)
        end
    end
end
            


end


colnames = tables_schema.names
coltypes = tables_schema.types

