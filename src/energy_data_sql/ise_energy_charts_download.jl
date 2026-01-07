"""
    download_ise_power_data(json_dir, year)

    json_dir : directory where to store downloaded json file
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
    path  = joinpath(json_dir, @sprintf("%s_%s.json", get_, year_))
    open(path, "w") do out
        JSON3.pretty(out, bb)
    end
    @info @sprintf("saved to %s", path)
end

"""
    Download ise data and store as json in json_dir
"""
function download_ise_data(json_dir, start_year, end_year)
    # cbpf Cross Border Physical Flows
    gets = ["public_power", "total_power", "installed_power", "cbpf"]
    for get in gets
        for year in start_year:end_year
            download_ise_power_data(json_dir, year, get)
        end
    end
end
