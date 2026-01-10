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

function DetrendedPowerData(df)
    DetrendedPowerData(
        df[!,:dates],
        df[!,:uts],
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
        df[!,:Bio_trend])
end

function scale_by_installed_power(power_data::DataFrame, installed_power::DataFrame)
    IP_Woff  = installed_power[!,:Woff ] ./ mean(installed_power[!,:Woff])
    IP_Won   = installed_power[!,:Won  ] ./ mean(installed_power[!,:Won])
    IP_Solar = installed_power[!,:Solar] ./ mean(installed_power[!,:Solar])
    IP_Bio   = installed_power[!,:Bio  ] ./ mean(installed_power[!,:Bio])

    Woff_ip  = @. power_data[!,:Woff ] / IP_Woff
    Won_ip   = @. power_data[!,:Won  ] / IP_Won
    Solar_ip = @. power_data[!,:Solar] / IP_Solar
    Bio_ip   = @. power_data[!,:Bio  ] / IP_Bio

    Woff_ip, Won_ip, Solar_ip, Bio_ip, Woff_trend, Won_trend, Solar_trend, Bio_trend
end

"""
    Woff, Won, Solar are scaled such that Woff+Won+Solar+Bio = Load
    Bio assumed not to increase further
"""
function scale_wind_and_solar_to_load!(Load, Woff, Won, Solar, Bio)
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


function get_detrended_power_data(power_data, installed_power, par)
    Woff_sc, Won_sc, Solar_sc, Bio_sc, Woff_trend, Won_trend, Solar_trend, Bio_trend  = 
        scale_by_installed_power(power_data, installed_power)

    Load_de, Load_trend  = detrend_time_series(power_data.Load)

    Woff_sc, Won_sc, Solar_sc, scale = scale_wind_and_solar_to_load!(Load_de, Woff_sc, Won_sc, Solar_sc, Bio_sc)
    @info @sprintf("scale_factor = %f", scale)

    WWSB_sc = Won_sc .+ Woff_sc .+ Solar_sc .+ Bio_de

    Load_sum = sum(Load_de)
    WWSB_sum = sum(WWSB_de)
    @info @sprintf("Load_sum = %10.4e, WWSB_sum = %10.4e, sum_L-sum_P = %10.4e", Load_sum, WWSB_sum, Load_sum-WWSB_sum)

    DataFrame([:dates, :uts, :Load, :Woff, :Won, :Solar, :Bio, :WWSBPower,
                :Load_trend, :Woff_trend, :Won_trend, :Solar_trend, :Bio_trend] 
                .=> 
                [power_data[!,:dates], power_data[!,:uts], Load_de, Woff_sc, Won_sc, Solar_sc, Bio_sc, WWSB_sc,
                Load_trend, Woff_trend, Won_trend, Solar_trend, Bio_trend]
    )
end

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
