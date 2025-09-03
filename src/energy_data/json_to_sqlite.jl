using Dates
using TimeZones
using SQLite
using Tables
using JSON3
using Printf
using Tables
using DataFrames
using CSV
using StringEncodings
using OrderedCollections
using SimpleLog

function get_jsonlists(jsonroot, categories)
    jlist = readdir(jsonroot)
    jfilelist = OrderedDict()
    for cat in categories
        jfilelist[cat] = sort(filter(x -> lowercase(splitext(x)[2]) == ".json" && occursin(cat, x), jlist))
    end
    jfilelist
end

"""
    create_umsaetze_table(db, table_name)

    create SQLite table
"""
function create_table(db, table_name, colnames)
    n = length(colnames)
    schema = Tables.Schema(colnames, types)
    SQLite.createtable!(db, table_name, schema; temp=false, ifnotexists=true)
end
"""
    create_unique_index(db, table_name)

    create unique index of table
"""
function create_unique_index(db, table_name)
    cnames, types = get_db_colnames_coltypes()
    nn = join(cnames, ", ")
    cmd = @sprintf("CREATE UNIQUE INDEX unique_all_cols ON %s (%s);", table_name, nn)
    DBInterface.execute(db, cmd)
end

function get_colnames(jfile)
    jdat = load_json(jfile)

    key1 = if haskey(jdat, "unix_seconds")
        "unix_seconds"
    elseif haskey(jdat, "time")
        "time"
    end

    key2 = if haskey(jdat, "production_types")
        "production_types"
    elseif haskey(jdat, "countries")
        "countries"
    end

    names = []
    for d in jdat[key2]
        push!(names, d["name"])
    end
    key1, key2, names
end

function normalize_colname(colname)
    colname = replace(colname, " " => "_")
    colname = replace(colname, "/" => "")
    colname = replace(colname, "__" => "_")
    colname = replace(colname, "-" => "_")
    colname = replace(colname, "_(incl._self_consumption)" => "")
    Symbol(colname)
end

function get_union_of_col_names(jlist)
    cols = Set()
    key1, key2 = "", ""
    for jfile in jlist
        key1, key2, colnames =  get_colnames(jfile)
        syms = normalize_colname.(colnames)
        for n in colnames
            sym = normalize_colname(n)
            push!(cols, sym)
            #println(n, " ", sym)
        end
    end
    cols = sort(collect(cols))
    cold = OrderedDict{Symbol, Int}()
    cold[Symbol(key1)] = 1
    j = 2
    for c in cols
        if !haskey(cold, c)
            cold[c] = j
            j +=1 
        end
    end
    cold
end

function load_data(jfile, colsyms0)
    key1, key2, colnames =  get_colnames(jfile)

    jdat = load_json(jfile);

    nr = length(jdat[key1])
    nc = length(colsyms0)

    vdat = Vector{Any}(undef,nc)
    vdat[1] = Vector{Int64}
    for i in 2:nc
        vdat[i] = zeros(Float64, nr)
    end

    j = colsyms0[Symbol(key1)]
    if j != 1
        @warne j
    end
    if key1 == "time"
        vdat[1] = parse.(Int64, jdat[key1])
    else
        vdat[1] = jdat[key1]
    end

    dats = jdat[key2];
    for (i,d) in enumerate(dats)
        dsym = normalize_colname(d["name"])
        j = colsyms0[dsym]
        vdat[j] = @. ifelse(d["data"] === nothing, 0.0, d["data"])
    end

    df = DataFrame(vdat, [k for (k,v) in colsyms0])
    df
end

function load_json(jsonfile)
    path = joinpath(jsonroot, jsonfile)
    data = open(path, "r") do io
        JSON3.read(read(io, String))
    end
    data
end

function load_files(jsonfiles, colsyms0)
    colsyms = collect(keys(colsyms0))
    coltypes = [Float64 for i in eachindex(colsyms)]
    coltypes[1] = Int64

    named_tuple = (; zip(colsyms, type[] for type in coltypes )...)
    df0 = DataFrame(named_tuple)
    for jsonfile in jsonfiles
        df = load_data(jsonfile, colsyms0)
        df0 = vcat(df0, df)
    end
    df0
end

function create_sqlite_db!(db, jsonfiles, categories)
    cat = "cbpf"
    for cat in categories
        colsyms0 = get_union_of_col_names(jsonfiles[cat])
        tt = if haskey(colsyms0, :unix_seconds)
            "unix_seconds"
        elseif haskey(colsyms0, :time)
            "time"
        else
            @warne colsyms0
        end
        @infoe cat, tt
        df0 = load_files(jsonfiles[cat], colsyms0)
        colnames = names(df0)

        tt = if "unix_seconds" in colnames
            "unix_seconds"
        elseif "time" in colnames
            "time"
        end
        if tt === nothing
            @infoe colnames
        end
        @infoe tt
        @infoe names(df0)
        @infoe SQLite.columns(db, "public_power")
        
        SQLite.load!(df0, db, cat)
        SQLite.removeduplicates!(db, cat, [tt])
        if :Load in colnames
            SQLite.execute(db, "DELETE FROM $cat WHERE Load < 1.0;")
        end
    end
end

function create_db()    
    dataroot = joinpath(dirname(dirname(@__DIR__)), "data")
    jsonroot = joinpath(dataroot, "json_downloads")
    #jsonroot = dataroot

    categories = ["cbpf", "public_power", "total_power", "installed_power"]
    jsonfiles = get_jsonlists(jsonroot, categories)
    dbname = joinpath(dataroot, "ise_data2.db")
    #rm(dbname)
    db = SQLite.DB(dbname)
    create_sqlite_db!(db, jsonfiles, categories)

    close(db)
end

function get_data(jahr_von, jahr_bis)
    dataroot = joinpath(dirname(dirname(@__DIR__)), "data")
    dbname = joinpath(dataroot, "ise_data.db")
    db = SQLite.DB(dbname)

    d1 = DateTime(jahr_von, 1, 1, 0, 0, 0)
    d2 = DateTime(jahr_bis, 12, 31, 23, 59, 59)
    uts1 = datetime2unix.(d1)
    uts2 = datetime2unix.(d2)
            
    cols = ["Load",
        "Biomass",
        "Solar",
        "Wind_offshore",
        "Wind_onshore",
        "Geothermal",
        "Hydro_Run_of_River",
        "Hydro_pumped_storage",
        "Hydro_pumped_storage_consumption",
        "Hydro_water_reservoir"]

    c = join(cols, ", ")

    total_power = DBInterface.execute(db, "SELECT $c FROM total_power WHERE unix_seconds >= $uts1 AND unix_seconds < $uts2;")
    installed_power = DBInterface.execute(db, "SELECT * FROM installed_power WHERE time >= $jahr_von AND time <= $jahr_bis;")

    DataFrame(total_power), DataFrame(installed_power)
end

#jahr_von = 2017
#jahr_bis = 2018
#total_power, installed_power = get_data(jahr_von, jahr_bis)

