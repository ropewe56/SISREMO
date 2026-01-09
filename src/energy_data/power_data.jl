function uconversion_factor(eto, efrom)
    @info eto, efrom
    Float64(getproperty(uconvert(eto, efrom), :val))
end

function average_to_hour(X)
    n = length(X)
    @. (X[1:4:n-3] + X[2:4:n-2] + X[3:4:n-1] + X[4:4:n]) * 0.25
end

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

"""
    EnergyData(eunit)

    EnergyData constructor
    par = PowerParameter()
"""
function PowerData(date1, date2, par)
    pp = select_from_db(date1, date2, ["public_power"])["public_power"]

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

    power_data = PowerData(dates, uts, Load, Woff, Won, Solar, Bio, Nuclear, WWSBPower)
    power_data
end

struct InstalledPowerData
    dates   :: Vector{DateTime} # 1
    uts     :: Vector{Int64}    # 2
    Woff    :: Vector{Float64}  # 3
    Won     :: Vector{Float64}  # 4
    Solar   :: Vector{Float64}  # 5
    Bio     :: Vector{Float64}  # 6
    BatCap  :: Vector{Float64}  # 7
    BatPow  :: Vector{Float64}  # 8start_year
end

"""
    uts_pubpower = power_data.uts
    punit = par.punit
"""
function iterpolate_installed_power(installed_power_y::DataFrame, punit, uts_pubpower)
    GW_to_unit = uconversion_factor(u_GW, 1.0*punit)

    year = installed_power_y[!,:time]
    uts = [datetime2unix(DateTime(@sprintf("%d-07-01", t))) for t in year]

    values = Vector{Any}([uts_pubpower])
    colnames = names(installed_power_y)
    for k in colnames[2:end]
        intp = linear_interpolation(uts, installed_power_y[!,k], extrapolation_bc=Line())
        push!(values, intp(uts_pubpower) .* GW_to_unit)
    end
    df = DataFrame(colnames .=> values)
    df
end

function minimum_value(v; scale = 1.0e-10)
    ipmin, ipmax = extrema(v)
    limit = scale * ipmax
    if ipmin < limit
        return @. ifelse(v < ipmin, limit, v)
    end
    v
end

function InstalledPowerData(power_data, par)
    date1, date2 = DateTime("2016-01-01"), now()
    installed_power_y = select_from_db(date1, date2, ["installed_power"])["installed_power"]

    installed_power_h = iterpolate_installed_power(installed_power_y, par.punit, power_data.uts)

    Woff    = minimum_value(installed_power_h[!,:Wind_offshore] * par.Woff_scale )
    Won     = minimum_value(installed_power_h[!,:Wind_onshore]  * par.Won_scale )
    SolarAC = minimum_value(installed_power_h[!,:Solar_AC]      * par.Solar_scale )
    SolarDC = minimum_value(installed_power_h[!,:Solar_DC]      * par.Solar_scale )
    Bio     = minimum_value(installed_power_h[!,:Biomass]       * par.Bio_scale )
    BatCap  = minimum_value(installed_power_h[!,:Battery_storage_capacity] )
    BatPow  = minimum_value(installed_power_h[!,:Battery_storage_power] )

    Solar = SolarAC .+ SolarDC

    InstalledPowerData(power_data.dates, power_data.uts, Woff, Won, Solar, Bio, BatCap, BatPow)
end

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

"""
    Woff, Won, Solar are scaled such that Woff+Won+Solar+Bio = Load
    Bio assumed not to increase further
"""
function scale_wind_and_solar!(Load, Woff, Won, Solar, Bio)
    WWS  = Woff .+ Won .+ Solar

    mean_Load = mean(Load)
    mean_Bio  = mean(Bio)
    mean_WWS  = mean(WWS)
    # mean_L*op = mean_R * scale + mean_B
    scale = (mean_Load - mean_Bio) / mean_WWS

    Woff  = Woff  .* scale
    Won   = Won   .* scale
    Solar = Solar .* scale

    Woff, Won, Solar, scale
end

function detrend_scaled_by_installed_power(power_data, installed_power)
    IP_Woff  = @. installed_power.Woff  / mean(installed_power.Woff)
    IP_Won   = @. installed_power.Won   / mean(installed_power.Won)
    IP_Solar = @. installed_power.Solar / mean(installed_power.Solar)
    IP_Bio   = @. installed_power.Bio   / mean(installed_power.Bio)

    Woff_ip  = @. power_data.Woff  / IP_Woff
    Won_ip   = @. power_data.Won   / IP_Won
    Solar_ip = @. power_data.Solar / IP_Solar
    Bio_ip   = @. power_data.Bio   / IP_Bio

    s1 = mean(power_data.Woff) / mean(Woff_ip)
    s2 = mean(power_data.Won)  / mean(Won_ip)
    s3 = mean(power_data.Solar)/ mean(Solar_ip)
    s4 = mean(power_data.Bio)  / mean(Bio_ip)

    Woff_de , Woff_trend  = detrend_time_series(Woff_ip  .* s1)
    Won_de  , Won_trend   = detrend_time_series(Won_ip   .* s2)
    Solar_de, Solar_trend = detrend_time_series(Solar_ip .* s3)
    Bio_de  , Bio_trend   = detrend_time_series(Bio_ip   .* s4)

    Woff_de, Woff_trend, Won_de, Won_trend, Solar_de, Solar_trend, Bio_de, Bio_trend 
end

function detrend_renewables(power_data)
    Woff_de , Woff_trend  = detrend_time_series(power_data.Woff)
    Won_de  , Won_trend   = detrend_time_series(power_data.Won)
    Solar_de, Solar_trend = detrend_time_series(power_data.Solar)
    Bio_de  , Bio_trend   = detrend_time_series(power_data.Bio)
    Woff_de, Woff_trend, Won_de, Won_trend, Solar_de, Solar_trend, Bio_de, Bio_trend 
end

function DetrendedPowerData(par, power_data, installed_power)    
    Woff_de, Woff_trend, Won_de, Won_trend, Solar_de, Solar_trend, Bio_de, Bio_trend  = if par.scale_with_installed_power_p
        scale_by_installed_power(power_data, installed_power)
    else
        detrend_renewables(power_data)
    end
    Load_de, Load_trend  = detrend_time_series(power_data.Load)

    Woff_sc, Won_sc, Solar_sc, scale = scale_wind_and_solar!(Load_de, Woff_de, Won_de, Solar_de, Bio_de)
    @info @sprintf("scale_factor = %f", scale)

    WWSB = Won_sc .+ Woff_sc .+ Solar_sc .+ Bio_de
    WWSB_de, WWSB_trend = detrend_time_series(WWSB)

    Load_sum = sum(Load_de)
    WWSB_sum = sum(WWSB_de)
    @info @sprintf("Load_sum = %10.4e, WWSB_sum = %10.4e, sum_L-sum_P = %10.4e", Load_sum, WWSB_sum, Load_sum-WWSB_sum)

    DetrendedPowerData(power_data.dates, power_data.uts,
        Load_de, Woff_sc, Won_sc, Solar_sc, Bio_de, WWSB_de,
        Load_trend, Woff_trend, Won_trend, Solar_trend, Bio_trend, WWSB_trend)
end

using Arrow

function save_detrended_power_data(detrended_power::DetrendedPowerData, detrended_power_path)
    colnames = ["uts", "Load", "Woff", "Won", "Solar", "Bio", "WWSBPower", 
                "Load_trend", "Woff_trend", "Won_trend", "Solar_trend", "Bio_trend", "WWSB_trend"]
    values = [detrended_power.uts,
        detrended_power.Load,
        detrended_power.Woff,
        detrended_power.Won,
        detrended_power.Solar,
        detrended_power.Bio,
        detrended_power.WWSBPower,
        detrended_power.Load_trend,
        detrended_power.Woff_trend,
        detrended_power.Won_trend,
        detrended_power.Solar_trend,
        detrended_power.Bio_trend,
        detrended_power.WWSB_trend]

    df = DataFrame(colnames .=> values)
    Arrow.write(detrended_power_path, df)
end

function load_detrended_power_data(detrended_power_path)
    df = Arrow.Table(detrended_power_path) |> DataFrame
    dates = unix2datetime.(df[!,:uts])
    DetrendedPowerData(dates, df[!,:uts], 
        df[!,:Load],
        df[!,:Woff],
        df[!,:Won],
        df[!,:Solar],
        df[!,:Bio],
        df[!,:WWSBPower],
        df[!,:Load_trend],
        df[!,:Woff_trend],
        df[!,:Won_trend],
        df[!,:Solar_trend],
        df[!,:Bio_trend],
        df[!,:WWSB_trend])
end

struct AveragedPowerData
    dates    :: Vector{DateTime} # 1
    uts      :: Vector{Int64}    # 2
    Load     :: Vector{Float64}  # 3
    Woff     :: Vector{Float64}  # 4
    Won      :: Vector{Float64}  # 5
    Solar    :: Vector{Float64}  # 6
    Bio      :: Vector{Float64}  # 7
#    Nuclear  :: Vector{Float64}  # 8
    WWSBPower:: Vector{Float64}  # 9 # sum of Woff Won Solar Bio
end

function AveragedPowerData(power_data, averaging_hours, averaging_method)
    Load     , dates_av = averaging(power_data.Load     , power_data.dates, averaging_hours, method = averaging_method)
    Woff     , dates_av = averaging(power_data.Woff     , power_data.dates, averaging_hours, method = averaging_method)
    Won      , dates_av = averaging(power_data.Won      , power_data.dates, averaging_hours, method = averaging_method)
    Solar    , dates_av = averaging(power_data.Solar    , power_data.dates, averaging_hours, method = averaging_method)
    Bio      , dates_av = averaging(power_data.Bio      , power_data.dates, averaging_hours, method = averaging_method)
#    Nuclear  , dates_av = averaging(power_data.Nuclear  , power_data.dates, averaging_hours, method = averaging_method)
    WWSBPower, dates_av = averaging(power_data.WWSBPower, power_data.dates, averaging_hours, method = averaging_method)

    uts_av = [Dates.datetime2unix(x) for x in dates_av]
    AveragedPowerData(dates_av, uts_av, Load, Woff, Won, Solar, Bio, WWSBPower)
end
