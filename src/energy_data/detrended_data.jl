using Dates

struct DetrendedPowerData
    dates       :: Vector{DateTime} # 1
    uts         :: Vector{Int64}    # 2
    Load        :: Vector{Float64}  # 3
    Woff        :: Vector{Float64}  # 4
    Won         :: Vector{Float64}  # 5 
    Solar       :: Vector{Float64}  # 6 
    Bio         :: Vector{Float64}  # 7 
    WWSBPower   :: Vector{Float64}  # 8 
    Load_trend  :: Vector{Float64}  # 9
    Woff_trend  :: Vector{Float64}  # 10 
    Won_trend   :: Vector{Float64}  # 11 
    Solar_trend :: Vector{Float64}  # 12 
    Bio_trend   :: Vector{Float64}  # 13 
    WWSB_trend  :: Vector{Float64}  # 14 
end

function save_detrended_power_data(dpd::DetrendedPowerData, hdf5_path)
    groups = Dict("detrended_power_data" => Dict(
        "uts"         => dpd.uts        ,
        "Load"        => dpd.Load       ,
        "Woff"        => dpd.Woff       ,
        "Won"         => dpd.Won        ,
        "Solar"       => dpd.Solar      ,
        "Bio"         => dpd.Bio        ,
        "WWSBPower"   => dpd.WWSBPower  ,
        "Load_trend"  => dpd.Load_trend ,
        "Woff_trend"  => dpd.Woff_trend ,
        "Won_trend"   => dpd.Won_trend  ,
        "Solar_trend" => dpd.Solar_trend,
        "Bio_trend"   => dpd.Bio_trend  ,
        "WWSB_trend"  => dpd.WWSB_trend ,        
    ))
    save_groups_as_hdf5(hdf5_path, groups, script_dir=false, permute_dims_p = true)
end

function load_detrended_power_data(hdf5_path)
    groups = load_groups_as_hdf5(hdf5_path, script_dir=false, permute_dims_p = true)
    dpd = groups["detrended_power_data"]

    uts         = dpd["uts"        ]
    Load        = dpd["Load"       ]
    Woff        = dpd["Woff"       ]
    Won         = dpd["Won"        ]
    Solar       = dpd["Solar"      ]
    Bio         = dpd["Bio"        ]
    WWSBPower   = dpd["WWSBPower"  ]
    Load_trend  = dpd["Load_trend" ]
    Woff_trend  = dpd["Woff_trend" ]
    Won_trend   = dpd["Won_trend"  ]
    Solar_trend = dpd["Solar_trend"]
    Bio_trend   = dpd["Bio_trend"  ]
    WWSB_trend  = dpd["WWSB_trend" ]
    
    dates = unix2datetime.(uts)

    DetrendedPowerData(dates, uts, Load, Woff, Won, Solar, Bio, WWSBPower, 
        Load_trend, Woff_trend, Won_trend, Solar_trend, Bio_trend ,WWSB_trend)
end

function save_detrended_power_data_to_db(dpd::DetrendedPowerData, dbname)
    DataFrame()
    colnames = [
        "uts"         ,
        "Load"        ,
        "Woff"        ,
        "Won"         ,
        "Solar"       ,
        "Bio"         ,
        "WWSBPower"   ,
        "Load_trend"  ,
        "Woff_trend"  ,
        "Won_trend"   ,
        "Solar_trend" ,
        "Bio_trend"   ,
        "WWSB_trend"  ] 

    values = [
        dpd.uts        ,
        dpd.Load       ,
        dpd.Woff       ,
        dpd.Won        ,
        dpd.Solar      ,
        dpd.Bio        ,
        dpd.WWSBPower  ,
        dpd.Load_trend ,
        dpd.Woff_trend ,
        dpd.Won_trend  ,
        dpd.Solar_trend,
        dpd.Bio_trend  ,
        dpd.WWSB_trend ]

    df = DataFrame(colnames .=> values)
    load_into_db(df, dbname, "detrended_power")    
end

function load_detrended_power_data_from_db(dbname)
    df = load_from_db(dbname, "detrended_power")

    uts         = df[!,:uts        ]
    Load        = df[!,:Load       ]
    Woff        = df[!,:Woff       ]
    Won         = df[!,:Won        ]
    Solar       = df[!,:Solar      ]
    Bio         = df[!,:Bio        ]
    WWSBPower   = df[!,:WWSBPower  ]
    Load_trend  = df[!,:Load_trend ]
    Woff_trend  = df[!,:Woff_trend ]
    Won_trend   = df[!,:Won_trend  ]
    Solar_trend = df[!,:Solar_trend]
    Bio_trend   = df[!,:Bio_trend  ]
    WWSB_trend  = df[!,:WWSB_trend ]
    
    dates = unix2datetime.(uts)

    DetrendedPowerData(dates, uts, Load, Woff, Won, Solar, Bio, WWSBPower, 
        Load_trend, Woff_trend, Won_trend, Solar_trend, Bio_trend, WWSB_trend)
end
