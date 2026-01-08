function create_table_power!(db, tablename)
    table_name = tablename
    unique_index = "unique_"*tablename

    SQLite.drop!(db, table_name; ifexists=true)
    SQLite.dropindex!(db, unique_index; ifexists=true)

    cmd = "CREATE TABLE $table_name
    (unix_seconds INT, 
    Load ,
    Residual_load ,
    Cross_border_electricity_trading ,
    Biomass ,
    Solar ,
    Wind_onshore ,
    Wind_offshore ,
    Hydro_Run_of_River ,
    Hydro_pumped_storage_consumption ,
    Hydro_pumped_storage ,
    Hydro_water_reservoir ,
    Waste ,
    Renewable_share_of_load ,
    Renewable_share_of_generation ,
    Fossil_gas ,
    Fossil_oil ,
    Geothermal ,
    Fossil_hard_coal ,
    Fossil_coal_derived_gas ,
    Fossil_brown_coal_lignite ,
    Nuclear ,
    Others );"

    DBInterface.execute(db, cmd)

    cmd = "CREATE UNIQUE INDEX $unique_index ON $table_name (unix_seconds);"
    DBInterface.execute(db, cmd)
end

function create_table_cbpf!(db)
    table_name = "cbpf"
    unique_index = "unique_cbpf"

    SQLite.drop!(db, table_name; ifexists=true)
    SQLite.dropindex!(db, unique_index; ifexists=true)

    cmd = "CREATE TABLE cbpf (
    unix_seconds INT NOT NULL, 
    Austria REAL,
    Belgium REAL,
    Czech_Republic REAL,
    Denmark REAL,
    France REAL,
    Luxembourg REAL,
    Netherlands REAL,
    Norway REAL,
    Poland REAL,
    Sweden REAL,
    Switzerland REAL);"

    DBInterface.execute(db, cmd)

    cmd = "CREATE UNIQUE INDEX $unique_index ON $table_name (unix_seconds);"
    DBInterface.execute(db, cmd)
end

function create_table_installed_power!(db)
    table_name = "installed_power"
    unique_index = "unique_installed_power"

    SQLite.drop!(db, table_name; ifexists=true)
    SQLite.dropindex!(db, unique_index; ifexists=true)

    cmd = "CREATE TABLE installed_power (
    time INT NOT NULL,
    Battery_storage_capacity REAL,
    Battery_storage_power REAL,
    Biomass REAL,
    Fossil_brown_coal_lignite REAL,
    Fossil_gas REAL,
    Fossil_hard_coal REAL,
    Fossil_oil REAL,
    Hydro REAL,
    Hydro_pumped_storage REAL,
    Nuclear REAL,
    Other_non_renewable REAL,
    Solar_AC REAL,
    Solar_DC REAL,
    Wind_offshore REAL,
    Wind_onshore REAL);"

    DBInterface.execute(db, cmd)

    cmd = "CREATE UNIQUE INDEX $unique_index ON $table_name (time);"
    DBInterface.execute(db, cmd)
end

function insert_data!(db, table_name, df)
    colnames = names(df)
    colnames = join(colnames, ",")

    DBInterface.execute(db, "BEGIN TRANSACTION;")
    try
        for row in eachrow(df)
            dat = collect(row)
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

function copy_from_db_to_db()
    dbname = joinpath(dataroot, "ise_data.sqlite")
    db = SQLite.DB(dbname)

    #create_table!(db, "total_power")
    #create_table!(db, "public_power")
    #create_table_cbpf!(db)
    #create_table_installed_power!(db)
    #SQLite.drop!(db, "installed_power"; ifexists=true)

    dbname2 = joinpath(dataroot, "ise_data.db")
    db2 = SQLite.DB(dbname2)

    df = DBInterface.execute(db2, "SELECT * FROM cbpf;") |> DataFrame
    df = select(df, Not(:sum))
    insert_data!(db, "cbpf", df)

    df = DBInterface.execute(db2, "SELECT * FROM total_power;") |> DataFrame
    insert_data!(db, "total_power", df)

    df = DBInterface.execute(db2, "SELECT * FROM public_power;") |> DataFrame
    insert_data!(db, "public_power", df)

    df = DBInterface.execute(db2, "SELECT * FROM installed_power;") |> DataFrame

    df = instdf["installed_power"]
    df_ = select(df, Not(:Wind_offshore_planned_WindSeeG, :Wind_onshore_planned_EEG_2023, :Solar_planned_EEG_2023))
    insert_data!(db, "installed_power", df_)
end