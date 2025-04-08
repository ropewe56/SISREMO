using Dates

struct AveragedPowerData
    dates    :: Vector{DateTime} # 1
    uts      :: Vector{Int64}    # 2
    Load     :: Vector{Float64}  # 3
    Woff     :: Vector{Float64}  # 4
    Won      :: Vector{Float64}  # 5
    Solar    :: Vector{Float64}  # 6
    Bio      :: Vector{Float64}  # 7
    Nuclear  :: Vector{Float64}  # 8
    WWSBPower:: Vector{Float64}  # 9
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