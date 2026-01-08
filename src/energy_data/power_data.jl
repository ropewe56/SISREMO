Base.@kwdef mutable struct PowerParameter
    punit                        = u_GW
    scale_Bio                    = 1.0
    Load_scale                   = 1.0
    Woff_scale                   = 1.0
    Won_scale                    = 1.0
    Solar_scale                  = 1.0
    Bio_scale                    = 1.0
    Geo_scale                    = 1.0
    SF1_factor                   = 1.0
    scale_to_installed_power_p   = true
    scale_with_installed_power_p = false
    averaging_hours              = 24*7
    averaging_method             = [:moving_average, :mean][1]
    plot_p                       = true
    plot_all_p                   = false
    log_p                        = false
    sisremo_root                 = SISREMOROOT
    data_dir                     = DATAROOT
    fig_dir                      = FIGDIR
    start_year                   = 2016
    end_year                     = 2025
end

struct PowerData
    dates     :: Vector{DateTime} # 1
    uts       :: Vector{Int64}    # 2
    Load      :: Vector{Float64}  # 3
    Woff      :: Vector{Float64}  # 4
    Won       :: Vector{Float64}  # 5
    Solar     :: Vector{Float64}  # 6
    Bio       :: Vector{Float64}  # 7
    Nuclear   :: Vector{Float64}  # 8
    WWSBPower :: Vector{Float64}  # 8
end

struct InstalledPowerData
    time   :: Vector{DateTime}  # 1
    Woff    :: Vector{Float64}  # 3
    Won     :: Vector{Float64}  # 4
    Solar   :: Vector{Float64}  # 5
    Bio     :: Vector{Float64}  # 6
    BatCap  :: Vector{Float64}  # 7
    BatPow  :: Vector{Float64}  # 8start_year
end

function InstalledPowerData(df)
    time   = df[!,:time]
    BatCap = df[!,:Battery_storage_capacity]
    BatPow = df[!,:Battery_storage_power]
    Bio    = df[!,:Biomass]
    Solar  = df[!,:Solar_AC] .+ df[!,:Solar_DC] 
    Woff   = df[!,:Wind_offshore]
    Won    = df[!,:Wind_onshore]
    InstalledPowerData(time, Woff, Won, Solar, Bio, BatCap, BatPow)
end


function uconversion_factor(eto, efrom)
    @info eto, efrom
    Float64(getproperty(uconvert(eto, efrom), :val))
end

function average_to_hour(X)
    n = length(X)
    @. (X[1:4:n-3] + X[2:4:n-2] + X[3:4:n-1] + X[4:4:n]) * 0.25
end

"""
    EnergyData(eunit)

    EnergyData constructor
    par = PowerParameter()
"""
function PowerData(date1, date2, par)
    data = select_from_db(date1, date2)
    pp = data["public_power"]

    # energy charts data are in MW, MW_to_unit is conversion factor to eunit (MW, GW, TW)
    MW_to_unit = uconversion_factor(par.punit, 1u_MW)
    
    Load    = average_to_hour(pp[!,:Load])          .* MW_to_unit .* par.Load_scale
    Woff    = average_to_hour(pp[!,:Wind_offshore]) .* MW_to_unit .* par.Woff_scale
    Won     = average_to_hour(pp[!,:Wind_onshore])  .* MW_to_unit .* par.Won_scale
    Solar   = average_to_hour(pp[!,:Solar])         .* MW_to_unit .* par.Solar_scale
    Bio     = average_to_hour(pp[!,:Biomass])       .* MW_to_unit .* par.Bio_scale
    Nuclear = average_to_hour(pp[!,:Nuclear])       .* MW_to_unit

    WWSBPower = Woff .+ Won .+ Solar .+ Bio

    uts = average_to_hour(pp[!,:unix_seconds])
    dates = Dates.unix2datetime.(uts)

    @info @sprintf("# timesteps = %d, length(Load) = %d, energy conversion = %e", length(uts), length(Load), MW_to_unit)

    df = data["installed_power"]
    names(df)

    PowerData(dates, uts, Load, Woff, Won, Solar, Bio, Nuclear, WWSBPower)
end


function load_ise_installed_power(hdf5_path)
    uts, prodtypes = load_installed_power_from_hdf5(hdf5_path)
    uts, prodtypes
end

function load_and_iterpolate_installed_power(hdf5_dir, punit, uts_pubpower)
    hp = @sprintf("installed_power.hdf5")
    hdf5_path = joinpath(hdf5_dir, hp)
    uts, installed_power = load_ise_installed_power(hdf5_path)

    GW_to_unit = uconversion_factor(u_GW, 1.0*punit)
    installed_power2 = Dict()

    for (i,(k,p)) in enumerate(installed_power)
        interp_linear_extrap = linear_interpolation(uts, p, extrapolation_bc=Line())
        installed_power2[k]  = interp_linear_extrap(uts_pubpower) .* GW_to_unit
    end
    @info "n_uts =", length(uts_pubpower)
    @info "size(installed_power2) =", length(keys(installed_power2))

    uts_pubpower, installed_power2
end

function InstalledPowerData(hdf5_dir, power_data, par)
    uts_pubpower, installed_power = load_and_iterpolate_installed_power(hdf5_dir, par.punit, power_data.uts)

    name_map = Dict("Bio"    => "Biomass",
                    "Woff"   => "Wind onshore",
                    "Won"    => "Wind offshore",
                    "Solar"  => "Solar",
                    "BatPow" => "Battery Storage (Power)",
                    "BatCap" => "Battery Storage (Capacity)")
    
    Woff    = installed_power[name_map["Woff"  ]] * par.Woff_scale
    Won     = installed_power[name_map["Won"   ]] * par.Won_scale
    Solar   = installed_power[name_map["Solar" ]] * par.Solar_scale
    Bio     = installed_power[name_map["Bio"   ]] * par.Bio_scale
    BatCap  = installed_power[name_map["BatCap"]]
    BatPow  = installed_power[name_map["BatPow"]]

    A = [Woff, Won, Solar, Bio, BatCap, BatPow]
    B = []
    n_rescaled = 0
    for a in A
        ipmin = minimum(a)
        ipmax = maximum(a)
        if ipmin < 1.0e-10*ipmax
            ipmin = 1.0e-10*ipmax
            a2 = @. ifelse(a < ipmin, ipmin, a)
            push!(B,  a2)
            n_rescaled += 1
        else
            push!(B, a)
        end
    end
    @info @sprintf("Installed power: # ipmin < 1.0e-10*ipmax = %d", n_rescaled)
    @info @sprintf("interpolated installed data %d", length(B[1]))
    #                                                    Woff  Won   Sol  Bio   BatCap BatPow
    InstalledPowerData(power_data.dates, power_data.uts, B[1], B[2], B[3], B[4], B[5], B[6])
end

struct AveragedPowerData
    dates    :: Vector{DateTime} # 1
    uts      :: Vector{Int64}    # 2
    Load     :: Vector{Float64}  # 3
    Woff     :: Vector{Float64}  # 4
    Won      :: Vector{Float64}  # 5
    Solar    :: Vector{Float64}  # 6
    Bio      :: Vector{Float64}  # 7
    Nuclear  :: Vector{Float64}  # 8
    WWSBPower:: Vector{Float64}  # 9 # sum of Woff Won Solar Bio
end

function get_averaged_power_data(power_data, averaging_hours, averaging_method)
    Load     , dates_av = averaging(power_data.Load     , power_data.dates, averaging_hours, method = averaging_method)
    Woff     , dates_av = averaging(power_data.Woff     , power_data.dates, averaging_hours, method = averaging_method)
    Won      , dates_av = averaging(power_data.Won      , power_data.dates, averaging_hours, method = averaging_method)
    Solar    , dates_av = averaging(power_data.Solar    , power_data.dates, averaging_hours, method = averaging_method)
    Bio      , dates_av = averaging(power_data.Bio      , power_data.dates, averaging_hours, method = averaging_method)
    Nuclear  , dates_av = averaging(power_data.Nuclear  , power_data.dates, averaging_hours, method = averaging_method)
    WWSBPower, dates_av = averaging(power_data.WWSBPower, power_data.dates, averaging_hours, method = averaging_method)

    uts_av = [Dates.datetime2unix(x) for x in dates_av]
    AveragedPowerData(dates_av, uts_av, Load, Woff, Won, Solar, Bio, Nuclear, WWSBPower)
end

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
