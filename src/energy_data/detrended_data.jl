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
