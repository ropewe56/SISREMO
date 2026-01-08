"""
    download_ise_power_data(json_dir, year)

    json_dir : directory where to store downloaded json file
    year : power data as json is downloaded for year and saved to power_<year>.json
"""
function download_ise_power_data(json_dir, year_, get_)
    query = [("country" => "de"), ("start", @sprintf("%s-01-01T00:00Z",year_)), ("end", @sprintf("%s-12-31T23:59Z", year_))]
    url   = "https://api.energy-charts.info/"*get_
    
    @info @sprintf("Download %s: %d", get_, year_)
    resp  = HTTP.get(url, ["Accept" => "application/json"], query = query)

    b     = String(resp.body)
    bb    = JSON3.read(b)
    path  = joinpath(json_dir, @sprintf("%s_%s.json", get_, year_))
    open(path, "w") do out
        JSON3.pretty(out, bb)
    end
    @info @sprintf("saved to %s", path)
end

"""
    download_ise_data(json_dir, start_year, end_year)
    
    Download ise data and store as json in json_dir

    public_power    Public Power
    total_power     Total Power
    cbpf            Cross Border Physical Flows
    installed_power Installed Power
"""
function download_ise_data(json_dir, start_year, end_year)
    gets = ["public_power", "total_power", "installed_power", "cbpf"]
    for get in gets
        for year in start_year:end_year
            download_ise_power_data(json_dir, year, get)
        end
    end
end

@inline function normalize_colname(colname)
    colname = replace(colname, " " => "_")
    colname = replace(colname, "/" => "")
    colname = replace(colname, "__" => "_")
    colname = replace(colname, "-" => "_")
    colname = replace(colname, "_(incl._self_consumption)" => "")
    colname = replace(colname, "(" => "")
    colname = replace(colname, ")" => "")
    colname = replace(colname, "," => "")
    colname
end

function normalize_keys(din)
    dout = OrderedDict()
    pk = collect(keys(din))
    for k in pk
        kn = normalize_colname(k)
        dout[kn] = din[k]
    end
    dout
end

"""
    json_to_dataframe(jsonpath)

read Energy Chart json file and convert it to a DataFrame

json keys are normalized for use DataFrames and SQLite

    total_power, public_power       : :unix_seconds, Vector of :production_types [:name, :data]
    cbpf (cross boundary power flow): :unix_seconds, Vector of :countries [:name, :data]
    installed_power                 : :time, Vector of :production_types [:name, :data]
"""
function json_to_dataframe(jsonpath)

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
    din = OrderedDict()
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
    DataFrame(colnames .=> values)
end

 