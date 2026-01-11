
function scale_detrend_print_info(data_ec, data)

    # dates
    dates = data_ec.dates
    @info @sprintf("start: (%d, %d, %d), stop: (%d, %d, %d)",
        Dates.day(dates[1]),   Dates.month(dates[1]),   Dates.year(dates[1]),
        Dates.day(dates[end]), Dates.month(dates[end]), Dates.year(dates[end]))

    # energy-chart data => 15 min time resolution, 1h = 60 * 60 * 1000 ms = 3.6e6 ms
    @info @sprintf("PE = %8.2e %sh", powers_to_energy_per_year(dates, P_ec), punit)
    @info @sprintf("LE = %8.2e %sh", powers_to_energy_per_year(dates, L_ec), punit)


    @info @sprintf("PE_de = %8.2e %sh", powers_to_energy_per_year(dates, P_de  ), punit)
    @info @sprintf("LE_de = %8.2e %sh", powers_to_energy_per_year(dates, L_de), punit)

    dates, L_ec, L, P_ec, P, P_de, P_trend, ΔEL, L_de, L_trend
end
function load_data(data_dir, punit, start_year, end_year; scale_to_installed_power = true, plot_p=false)

    data = PowerData(data_dir, punit, start_year, end_year);
    L_ec = data.Load
    dates = data.dates

    # energy charts
    P_ec = @. data.Woff + data.Won + data.Solar + data.Bio;


    # dates
    @info @sprintf("start: (%d, %d, %d), stop: (%d, %d, %d)",
        Dates.day(dates[1]),   Dates.month(dates[1]),   Dates.year(dates[1]),
        Dates.day(dates[end]), Dates.month(dates[end]), Dates.year(dates[end]))

    # energy-chart data => 15 min time resolution, 1h = 60 * 60 * 1000 ms = 3.6e6 ms
    @info @sprintf("PE = %8.2e %sh", powers_to_energy_per_year(dates, P_sc), punit)
    @info @sprintf("LE = %8.2e %sh", powers_to_energy_per_year(dates, L_sc), punit)

    P_de, P_trend, ΔEL, L_de, L_trend = scale_and_detrend(L_ec, P_ec);

    @info @sprintf("PE_de = %8.2e %sh", powers_to_energy_per_year(dates, P_de), punit)
    @info @sprintf("LE_de = %8.2e %sh", powers_to_energy_per_year(dates, L_de), punit)


    dates, data.Load, L, P_ec, P, P_de, P_trend, ΔEL, L_de, L_trend
end

"""
    scale_and_detrend(Load, RP)

    renewable energy RP is scaled such that mean(RP) == mean(Load) over whole time considered (e.g. 2015 - 2022)
"""
function scale_and_detrend(Load::Vector{Float64}, RP::Vector{Float64})
    n = length(Load)
    nh = div(n,2)

    # determine trend by fitting data with a k'th order polynomial
    k = 2

    # detrend Load
    Load_trend = polynomial_fit(Load, k)
    Load_de = @. Load * Load_trend[nh] / Load_trend

    scale = (mean(Load) / mean(RP))
    #@info @sprintf("(mean(RP) / mean(Load)) = %3.2f", 1.0/scale)

    # scale RP
    RP_sc = @. RP * scale

    # detrend RP
    RP_trend = polynomial_fit(RP_sc, k)
    RP_de = @. RP_sc * RP_trend[nh] / RP_trend

    # scale RP_de using detrended data
    RP_de_sc = RP_de * mean(Load_de)/mean(RP_de)

    # diff between detrended RP and detrended Load
    ΔEL = (RP_de_sc - Load_de)

    RP_de_sc, RP_trend, ΔEL, Load_de, Load_trend
end
