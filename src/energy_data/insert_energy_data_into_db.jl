function get_last_unix_seconds(table, db)
    df = DBInterface.execute(db, "SELECT unix_seconds FROM $table;") |> DataFrame
    maximum(df[!,:unix_seconds])
end

function filter_production_after(df, table, db)
    last_uts = get_last_unix_seconds(table, db)
    filter(row -> row[:unix_seconds] > last_uts, df)
end

function downloaded_json_to_dataframes(tables, year)
    alldf = Dict()
    for table in tables
        jsonpath = joinpath(jsonroot, @sprintf("%s_%s.json", table, year))
        df = json_to_dataframe(jsonpath)
        alldf[table] = df
    end
    alldf
end

function add_new_data_to_db()
    dbname = joinpath(dataroot, "ise_data.sqlite")
    db = SQLite.DB(dbname)

    tables = ["installed_power"]
    instpower = downloaded_json_to_dataframes(tables, 2025)

    tables = ["cbpf", "public_power", "total_power"]
    allpower = downloaded_json_to_dataframes(tables, 2025)

    table = "cbpf"
    df = select(allpower[table], Not(:sum))
    insert_data!(db, table, df)
    
    table = "public_power"
    df = allpower[table]
    insert_data!(db, table, df)

    table = "total_power"
    df = allpower[table]
    insert_data!(db, table, df)

    df = instpower["installed_power"]
    df_ = select(df, Not(:Wind_offshore_planned_WindSeeG, :Wind_onshore_planned_EEG_2023, :Solar_planned_EEG_2023))
    insert_data!(db, "installed_power", df_)
end

#start_year, end_year = dataroot, 2025, 2025
#download_ise_data(jsonroot, start_year, end_year)
#add_new_data_to_db()

