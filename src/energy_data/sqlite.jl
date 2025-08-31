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

function get_json_list(root, prefix)
    jsonfiles = []
    for (rr,dd,ff) in walkdir(root)
        csv = filter(x -> lowercase(splitext(x)[2]) == ".json" && occursin(prefix, x), ff)
        for c in csv
            push!(jsonfiles, joinpath(rr, c))
        end
    end
    jsonfiles
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
    Symbol(colname)
end

function get_union_of_col_names(jlist)
    cols = Set()
    for jfile in jlist
        key1, key2, colnames =  get_colnames(jfile)
        syms = normalize_colname.(colnames)
        for n in syms
            push!(cols, n)
        end
    end
    cols = sort(collect(cols))
    cold = OrderedDict{Symbol, Int}()
    cold[:unix_seconds] = 1
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
    vdat = Vector{Vector{Float64}}(undef,nc)
    for i in 1:nc
        vdat[i] = zeros(Float64, nr)
    end

    j = colsyms0[Symbol(key1)]
    vdat[j] = jdat[key1]

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
    named_tuple = (; zip(colsyms, type[] for type in coltypes )...)

    df0 = DataFrame(named_tuple)
    for jsonfile in jsonfiles
        df = load_data(jsonfile, colsyms0)
        df0 = vcat(df0, df)
    end
    df0
end


colsyms0 = get_union_of_col_names(jsonfiles["public_power"])
df0 = load_files(jsonfiles["public_power"], colsyms0)
db = SQLite.DB("public_power.db")
SQLite.load!(df0, db, "public_power")
SQLite.removeduplicates!(db, "public_power", ["unix_seconds"])
SQLite.execute(db, "DELETE FROM public_power WHERE Load == 0.0;")


