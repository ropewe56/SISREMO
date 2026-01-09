function select_from_db(date1, date2, tables)
    db = SQLite.DB(DBPATH)

    uts1 = datetime2unix(date1)
    uts2 = datetime2unix(date2)

    year1 = Dates.year(date1)
    year2 = Dates.year(date2)

    d = Dict()
    for table in tables
        if table == "installed_power"
            cmd = "SELECT * FROM installed_power WHERE time >= $year1 AND time <= $year2 ORDER BY time;"
            d["installed_power"] = DBInterface.execute(db, cmd) |> DataFrame
        else
            cmd = "SELECT * FROM $table WHERE unix_seconds >= $uts1 AND unix_seconds <= $uts2 ORDER BY unix_seconds;"
            d[table] = DBInterface.execute(db, cmd) |> DataFrame
        end
    end
    d
end
