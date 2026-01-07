include("include_energy_data.jl")

function get_last_unix_seconds(table, db)
    df = DBInterface.execute(db, "SELECT unix_seconds FROM $table;") |> DataFrame
    maximum(df[!,:unix_seconds])
end

function filter_production_after(alldf, db)
    alldf_after = OrderedDict()
    for tbl in ["public_power", "total_power", "cbpf"]
        last_uts = get_last_unix_seconds(tbl, db)
        alldf_after[tbl] = filter(row -> row[:unix_seconds] > last_uts, alldf[tbl])
    end
    alldf_after
end

function get_tables(df, db)
    function insert_entry!(d, k, e)
        if haskey(d, k)
            push!(d[k], e)
        else
            d[k] = [e]
        end
    end

    tblschema = OrderedDict()
    for table in SQLite.tables(db)
        if haskey(df, table.name)
            dbcolnames = Symbol.(table.schema.names)
            dfcolnames = Symbol.(names(df[table.name]))
            for cn in dfcolnames
                if cn in dbcolnames
                    insert_entry!(tblschema, table.name, cn)
                end
            end
            for cn in dbcolnames
                if !(cn in dfcolnames)
                    @warn cn, "not in", dfcolnames
                end
            end
        end
    end
    tblschema
end

#start_year, end_year = dataroot, 2025, 2025
#download_ise_data(jsonroot, start_year, end_year)

tables = ["cbpf", "public_power", "total_power", "installed_power"]
function json_to_dataframes(tables, year)
    alldf = Dict()
    for table in tables
        jsonpath = joinpath(jsonroot, @sprintf("%s_%s.json", table, year))
        df = json_to_dataframe(jsonpath)
        alldf[table] = df
    end
    alldf
end

proddf = json_to_dataframes(["cbpf", "public_power", "total_power"], 2025)
instdf = Dict("installed_power" => json_to_dataframe(joinpath(jsonroot, @sprintf("%s_%d.json", "installed_power", 2025))))


dbname = joinpath(dataroot, "ise_data.db")
db = SQLite.DB(dbname)

df_newest = filter_production_after(proddf, db)

function insert_data!(db, table_name, df)
    colnames = names(df)
    push!(colnames, "Nuclear")
    colnames = join(colnames, ",")

    DBInterface.execute(db, "BEGIN TRANSACTION;")
    try
        for row in eachrow(df)
            dat = collect(row)
            push!(dat, 0.0)
            fz = join(repeat("?", length(dat)), ",")
            #@info "INSERT INTO $table_name ($colnames) VALUES ($fz)"
            SQLite.execute(db, "INSERT INTO $table_name ($colnames) VALUES ($fz)", dat)
        end
        DBInterface.execute(db, "COMMIT TRANSACTION;")
    catch e
        @warn e
        DBInterface.execute(db, "ROLLBACK TRANSACTION")
    end
end

for (table_name, df) in df_newest
    table_name = "total_power"
    df = df_newest[table_name]
    insert_data!(db, table_name, df)
end

unix2datetime(get_last_unix_seconds("public_power", db))
unix2datetime(get_last_unix_seconds("total_power", db))
unix2datetime(get_last_unix_seconds("cbpf", db))

df_newest



names = df_newest["public_power"]
for row in eachrow(df_newest["public_power"])
    "INSERT INTO public_power (names) VALUES()"
end

tblschema = get_tables(alldf_after, db)

for tb in tblschema["public_power"]
    println(tb)
end


for n in names(alldf_after["public_power"])
    println(n)
end


colnames = tables_schema.names
coltypes = tables_schema.types

columns = [
:unix_seconds,
:Load,
:Residual_load,
:Cross_border_electricity_trading,
:Biomass,
:Solar,
:Wind_onshore,
:Wind_offshore,
:Hydro_Run_of_River,
:Hydro_pumped_storage_consumption,
:Hydro_pumped_storage,
:Hydro_water_reservoir,
:Waste,
:Renewable_share_of_load,
:Renewable_share_of_generation,
:Fossil_gas,
:Fossil_oil,
:Geothermal,
:Fossil_hard_coal,
:Fossil_coal_derived_gas,
:Fossil_brown_coal_lignite,
:Nuclear,
:Others,
]
