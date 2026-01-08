function select_from_db(date1, date2)
    db = SQLite.DB(DBPATH)

    uts1 = datetime2unix(date1)
    uts2 = datetime2unix(date2)

    year1 = year(date1)
    year2 = year(date2)

    d = Dict()
    tables = ["total_power", "public_power", "cbpf"]
    for table in tables
        cmd = "SELECT * FROM $table WHERE unix_seconds >= $uts1 AND unix_seconds <= $uts2 ORDER BY unix_seconds;"
        d[table] = DBInterface.execute(db, cmd) |> DataFrame
    end
    cmd = "SELECT * FROM installed_power WHERE time >= $year1 AND time <= $year2 ORDER BY time;"
    d["installed_power"] = DBInterface.execute(db, cmd) |> DataFrame

    d
end

#date1 = DateTime("2020-01-01")
#date2 = DateTime("2022-12-31")
#d = select_from_db(date1, date2)