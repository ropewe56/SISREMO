using Common
using Statistics
using Interpolations

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("storage/utils.jl")
#include("storage/hdf5_utils.jl")
include("storage/data_energy.jl")

root = dirname(@__DIR__)
get_data_dir()   = joinpath(root, "data")
get_fig_dir()    = joinpath(root, "figures")
mkpath(get_fig_dir())

"""
    get the elapsed time of 1 step
    Dates.value(DateTime) => ms since 1 AD
    energy-chart data => 15 min time resolution ΔTh = 0.25
    1h = 60 * 60 * 1000 ms = 3.6e6 ms
    returns ΔTh - multiple of an hour per stpe
"""
function get_step_ΔTh(dates)
    ms = Dates.value(dates[2] - dates[1])
    ms/(3.6e6)
end

"""
    Get the number of years of the time series given by dates
"""
function number_years(dates)
    ΔT_ms_e  = Dates.value(dates[end] - dates[1])
    Δh_e     = ΔT_ms_e/(3.6e6)
    Δh_e/(365*24)
end

"""
    ehergy produced per year by power series P
"""
function powers_to_energy_per_year(dates, P)
    ΔTh = get_step_ΔTh(dates)
    nb_years = number_years(dates)
    sum(P)*ΔTh/nb_years
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
    #@infoe @sprintf("(mean(RP) / mean(Load)) = %3.2f", 1.0/scale)

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

mutable struct Storage
    SC   :: Float64         # storage capacity [punit_h] TWh
    SF   :: Vector{Float64} # storage fill     [punit_h] TWh
    ηin  :: Float64         # efficiency while charging
    ηout :: Float64         # efficiency while sdischarging
    I2   :: Vector{Float64} # direct flow to load           [punit] TW, GW, MW
    I3   :: Vector{Float64} # to storage and/or curtailment [punit] TW, GW, MW
    I4   :: Vector{Float64} # to storage                    [punit] TW, GW, MW
    I5   :: Vector{Float64} # to curtailment                [punit] TW, GW, MW
    I6   :: Vector{Float64} # from storage to load          [punit] TW, GW, MW
    I7   :: Vector{Float64} # other sources to load         [punit] TW, GW, MW
end
function Storage(SC::Float64, nb_steps::Int64; ηin = 1.0, ηout = 1.0)
    SF = zeros(Float64, nb_steps)
    I2 = zeros(Float64, nb_steps)
    I3 = zeros(Float64, nb_steps)
    I4 = zeros(Float64, nb_steps)
    I5 = zeros(Float64, nb_steps)
    I6 = zeros(Float64, nb_steps)
    I7 = zeros(Float64, nb_steps)
    Storage(SC, SF, ηin, ηout, I2, I3, I4, I5, I6, I7)
end

function write_to_log(stg::Storage, ΔPin, out, i, j)
    B = stg.IF[i] - stg.I6[i] - (stg.SF[i] - stg.SF[i-1])
    write(out, @sprintf("%5d  %d  ΔP = %9.2e, S = %9.2e, I = %9.2e, O = %9.2e, C = %9.2e, M = %9.2e, B = %9.2e\n", i, j, ΔPin, stg.SF[i], stg.IF[i], stg.OF[i], stg.CT[i], stg.RE[i], B))
end

function write_power_step_to_log(stg::Storage, out, i, j)
    write(out, @sprintf("%5d  %d  P = %9.2e, L = %9.2e, I2 = %9.2e, I3 = %9.2e, I4 = %9.2e, I5 = %9.2e, I6 = %9.2e, I7 = %9.2e\n",
        i, j, stg.SF[i], P, L, stg.I2[i], stg.I3[i], stg.I4[i], stg.I5[i], stg.I6[i], stg.I7[i]))
end

function power_step(stg::Storage, L, P, i; out = nothing)
    ΔSC = stg.SC - stg.SF[i-1]
    if P - L > 0.0
        I2[i] = L
        I3[i] = P - L
        I4[i] = min(I3[i], ΔSC/stg.ηin)
        stg.SF[i] = stg.SF[i-1] + I4[i]*stg.ηin
        I5[i] = max(0.0, I3[i]-I4[i])  # curtailment
    else
        I2[i] = P
        I6[i] = min(stg.SF[i-1] * stg.ηout, L-P)
        stg.SF[i] = stg.SF[i-1] - I6[i]/stg.ηout
        I7[i] = L - I2[i] - I6[i]      # residual load
    end
end

"""
    compute_storage_level(dates, Load, RP, punit, over_production, storage_capacity)

    given Load, RP, over_production and storage_capacity compute storage level as a funtion of time

    dates : times, ΔT = 1h
    Load  : power consumed
    RP    : renewable power production
    punit : power unit of Load and RP (MW, GW, TW)
    over_production  : renewable over production capacity factor, 1.0 is no over production capacity
    storage_capacity : storage capacity
"""
function compute_storage_level(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64},
    st_capacities, oprod::Float64; do_log = false)

    # renewable energy over production
    RP = RP .* oprod

    nb_steps  = size(RP, 1)
    for stc in st_capacities
        storage = Storage(stc, nb_steps)
        storage.SF[1] = stc*0.5
        push!(storages, storage)
    end

    for i in 2:n
        power_step(storages[1], L, P, i, out = out1)
        for j in 2:length(storages)
            charge(storages[j], storages[j-1].I7[i], storages[j-1].I5[i], i, out = out2)
        end
    end

    for (j, stg) in enumerate(storages)
        direct = sum(stg.I2)
        st_in  = sum(stg.I4)
        st_out = sum(stg.I6)
        other  = sum(stg.I7)
        curt   = sum(stg.I5)
        @infoe @sprintf("j = %d, op = %3.1f, SC = %8.2e, SF[end] = %8.2e", j, oprod, stg.SC, stg.SF[end])
        @infoe @sprintf("direct = %8.2e, st_in = %8.2e, st_out = %8.2e, other = %8.2e, curt = %8.2e", direct, st_in, st_out, other, curt)
    end

    storages
end

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function compute_storage_fill_level(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64}, punit::String, stc1, stc2, oprod; do_log = false)
    @infoe @sprintf("==== compute_storage_fill_level ===================================================================")

    res1 = []
    res2 = []
    for (sc1, sc2, op) in zip(stc1, stc2, oprod)
        stg1, stg2 = compute_storage_level(dates, Load, RP, sc1, sc2, op, do_log = do_log)
        push!(res1, stg1)
        push!(res2, stg2)
    end
    (res1, res2, oprod)
end

"""
    determine_overproduction(Load, RP)

    determine minimum storage capacity as a function of over production

    Load : detrended Load
    RP : detrended and scaled renewable power
    punit : unit of Load and RP (MW. GW, TW)
"""
function determine_overproduction(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64}, punit::String)
    overproduction = collect(LinRange(1.05, 1.5, 20))
    storage_capacities = []
    for op in overproduction
        stc1 = 1.0
        stc2 = 1.0
        minS = -1.0
        it = 0
        while minS < 0.0 && it < 50
            stg1, stg2 = compute_storage_level(dates, Load, RP, op, stc1, stc2)

            min_storage_level1 = minimum(stg1.SF)
            min_storage_level2 = minimum(stg2.SF)

            stc1 = stc1 - min_storage_level1
            stc2 = stc2 - min_storage_level2
            it += 1
        end
        push!(storage_capacities, (stc1, stc2))
    end
    pl.plot(overproduction, storage_capacities, "r.")
end

"""
    plot_powers(dates, Load, RP, averaging_hours, fig_dir, punit, fig)

    plot Load and RP

    dates - times
    Load -
    RP   -
    averaging_hours - if data are averaged number of hours to avergae over
    fig_dir - directory where figures are saved to
    punit - (MW, GW, TW)
    fig - number of matplotlib figure
"""
function plot_powers(dates::Vector{DateTime}, Load_ec::Vector{Float64}, Load::Vector{Float64}, RP_ec::Vector{Float64}, RP::Vector{Float64},
    averaging_hours::Int64, fig_dir::String, punit::String, fig::Vector{Int64})

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, Load_ec, label="Load")
    pl.plot(dates, RP_ec, label="Renewables")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", punit))
    pl.grid()
    pl.title("Energy-Charts data")
    pl.legend()

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, RP, label="Renewables")
    pl.plot(dates, Load, label="Load")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", punit))
    pl.grid()
    pl.title("scaled  data")
    pl.legend()

    if averaging_hours > 1
        pl.title(@sprintf("averaged, %d days window", div(averaging_hours,24)))
        pl.savefig(joinpath(fig_dir, @sprintf("RP_averaged.png")))
    else
        @infoe fig_dir
        pl.savefig(joinpath(fig_dir, @sprintf("RP.png")))
    end
end

function plot_detrended(dates::Vector{DateTime}, RP::Vector{Float64}, RP_de::Vector{Float64}, RP_trend::Vector{Float64}, ΔEL::Vector{Float64},
    Load::Vector{Float64}, Load_de::Vector{Float64}, Load_trend::Vector{Float64},
    punit ::String, fig_dir::String, fig::Vector{Int64}; data_are_averaged = false)

    label1 = "Load"
    label2 = "Load_de"
    path1  = "Load_detrended.png"
    path2  = "Load_diff_detrended.png"
    if data_are_averaged
        label1 = "Load_av"
        label2 = "Load_av_de"
        path1  = "Load_av_detrended.png"
    end

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, Load      , label = label1)
    pl.plot(dates, Load_de   , label = label2)
    pl.plot(dates, Load_trend, "r", linewidth=3, label = "Load_trend")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", punit))
    pl.grid()
    pl.legend()
    pl.title("Load detrended")
    pl.savefig(joinpath(fig_dir, path1))

    label1 = "RP"
    label2 = "RP_de"
    label3 = "RP_de - Load"
    path1  = "RP_detrended.png"
    path2  = "RP_diff_detrended.png"
    if data_are_averaged
        label1 = "RP_av"
        label2 = "RP_av_de"
        label3 = "RP_av_de - Load_av"
        path1  = "RP_av_detrended.png"
        path2  = "RP_av_diff_detrended.png"
    end

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, RP      , label = label1)
    pl.plot(dates, RP_de   , label = label2)
    pl.plot(dates, RP_trend,  "r", linewidth=3, label = "RP_trend")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", punit))
    pl.grid()
    pl.legend()
    pl.title("Renewable detrended")
    pl.savefig(joinpath(fig_dir, path1))

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, ΔEL, label = label3)
    pl.xlabel("time")
    pl.ylabel(@sprintf("(P_E - P_L) [%s]", punit))
    pl.legend()
    pl.grid()
    pl.title("Renewable power - Load")
    pl.savefig(joinpath(fig_dir, path2))
end

function plot_storage_fill_level(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64},
    storage, oprod, title, fig_dir::String, fig::Vector{Int64}, pngpath::String, punit ::String; plot_all_p = false)

    prozent = "%"
    pl.figure(fig[1]); fig[1] += 1
    for (stg, op) in zip(storage, oprod)
        label = @sprintf("%s storage_cpacity=%2.fh %s, op = %3.f %s", title, stg.SC, punit, (op-1.0)*100.0, prozent)
        pl.plot(dates, stg.SF, label = label)
    end
    pl.xlabel("time")
    pl.ylabel(@sprintf("storage fill level [%sh]", punit))
    pl.legend()
    pl.grid()
    pl.savefig(joinpath(fig_dir, pngpath.*".png"))

    if plot_all_p
        for (stg, op) in zip(storage, oprod)
            sc = stg.SC

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, RP .* op ,   label="RP")
            pl.plot(dates, Load,   label="Load")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Load, RP: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_1_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Load,   label="Load")
            pl.plot(dates, stg.IF,   label="IF")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Storage power inflow: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_2_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Load,   label="Load")
            pl.plot(dates, stg.OF,  label="OF")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Storage power outflow: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_3_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Load,   label="Load")
            pl.plot(dates, stg.IM, label="IM")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Power import: SC=%2.f[%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_4_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Load,   label="Load")
            pl.plot(dates, stg.CT, label="CT")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Power curtailment: SC=%2.f %s, op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_5_%2.f_%1.f", title, sc, (op-1.0)*100.0)))
        end
    end
end

function load_data(data_dir, punit, start_year, end_year; scale_to_installed_power = true, plot_p=false)
    data = PowerData(data_dir, punit, start_year, end_year);
    Load = data.Load

    local RP

    # energy charts
    RP_ec = @. data.Woff + data.Won + data.Solar + data.Bio;

    if scale_to_installed_power
        IP = InstalledPower(data, data_dir, punit, start_year, end_year)

        md_Woff  = 1.0 / mean(IP.Woff)
        md_Won   = 1.0 / mean(IP.Won)
        md_Solar = 1.0 / mean(IP.Solar)
        md_Bio   = 1.0 / mean(IP.Bio)

        IP_Woff  = @. IP.Woff  * md_Woff
        IP_Won   = @. IP.Won   * md_Won
        IP_Solar = @. IP.Solar * md_Solar
        IP_Bio   = @. IP.Bio   * md_Bio

        Woff  = @. data.Woff  / IP_Woff
        Won   = @. data.Won   / IP_Won
        Solar = @. data.Solar / IP_Solar
        Bio   = @. data.Bio   / IP_Bio

        # renewables is sum of wind onshore, wind offshore and solar
        RP0 = @. Woff + Won + Solar + Bio;

        scale = mean(Load) / mean(RP0)
        @infoe mean(Load), mean(RP0), scale
        RP = @. RP0 * scale

        if plot_p
            plt.figure()
            plt.plot(data.dates, IP.Solar, label="Solar")
            plt.plot(data.dates, IP.Won,   label="Won")
            plt.plot(data.dates, IP.Woff,  label="Woff")
            plt.plot(data.dates, IP.Bio,   label="Bi ")
            plt.legend()
            plt.title("IP 1")

            plt.figure()
            plt.plot(data.dates, IP_Solar, label="Solar")
            plt.plot(data.dates, IP_Won,   label="Won")
            plt.plot(data.dates, IP_Woff,  label="Woff")
            plt.plot(data.dates, IP_Bio,   label="Bi ")
            plt.legend()
            plt.title("IP 2")

            plt.figure()
            plt.plot(data.dates, data.Solar, label="Solar")
            plt.plot(data.dates, data.Won,   label="Won")
            plt.plot(data.dates, data.Woff,  label="Woff")
            plt.plot(data.dates, data.Bio,   label="Bi ")
            plt.legend()
            plt.title("EE 1")

            plt.figure()
            plt.plot(data.dates, Solar, label="Solar")
            plt.plot(data.dates, Won,   label="Won")
            plt.plot(data.dates, Woff,  label="Woff")
            plt.plot(data.dates, Bio,   label="Bi ")
            plt.legend()
            plt.title("EE 2")

            plt.figure()
            plt.plot(data.dates, RP0)
            plt.title("RP 1")

            plt.figure()
            plt.plot(data.dates, RP)
            plt.title("RP 2")

            plt.figure()
            plt.plot(data.dates, Load)
            plt.title("Load")
        end
    else
        RP = @. data.Woff + data.Won + data.Solar + data.Bio;
    end

    # dates
    dates = data.dates
    @infoe @sprintf("start: (%d, %d, %d), stop: (%d, %d, %d)",
        Dates.day(dates[1]),   Dates.month(dates[1]),   Dates.year(dates[1]),
        Dates.day(dates[end]), Dates.month(dates[end]), Dates.year(dates[end]))

    # energy-chart data => 15 min time resolution, 1h = 60 * 60 * 1000 ms = 3.6e6 ms
    @infoe @sprintf("RE = %8.2e %sh", powers_to_energy_per_year(dates, RP  ), punit)
    @infoe @sprintf("LE = %8.2e %sh", powers_to_energy_per_year(dates, Load), punit)

    RP_de, RP_trend, ΔEL, Load_de, Load_trend = scale_and_detrend(Load, RP);

    @infoe @sprintf("RE_de = %8.2e %sh", powers_to_energy_per_year(dates, RP_de  ), punit)
    @infoe @sprintf("LE_de = %8.2e %sh", powers_to_energy_per_year(dates, Load_de), punit)


    dates, data.Load, Load, RP_ec, RP, RP_de, RP_trend, ΔEL, Load_de, Load_trend
end

"""
    load data and compute and plot storage fille levels, original times (15 min)
"""
function comp_and_plot(stc1, stc2, oprod, data_dir, fig_dir, punit, start_year, stop_year; plot_p = false, plot_all_p = false, do_log = false)

    dates, Load_ec, Load,RP_ec, RP, RP_de, RP_trend, ΔEL, Load_de, Load_trend = load_data(data_dir, punit, start_year, stop_year)

    (res1, res2, oprod) = compute_storage_fill_level(dates, Load_de, RP_de, punit, stc1, stc2, oprod, do_log = do_log)

    for (stg1, stg2, op) in zip(res1, res2, oprod)
        Load_per_year = powers_to_energy_per_year(dates, Load)
        IF1  = powers_to_energy_per_year(dates, stg1.IF) / Load_per_year
        OF1  = powers_to_energy_per_year(dates, stg1.OF) / Load_per_year
        CT1  = powers_to_energy_per_year(dates, stg1.CT) / Load_per_year
        IM1  = powers_to_energy_per_year(dates, stg1.IM) / Load_per_year
        IF2  = powers_to_energy_per_year(dates, stg2.IF) / Load_per_year
        OF2  = powers_to_energy_per_year(dates, stg2.OF) / Load_per_year
        CT2  = powers_to_energy_per_year(dates, stg2.CT) / Load_per_year
        IM2  = powers_to_energy_per_year(dates, stg2.IM) / Load_per_year
        SF1  = stg1.SF[end] / Load_per_year
        SF2  = stg2.SF[end] / Load_per_year

        @infoe @sprintf("==== comp_and_plot ================================================================================")
        @infoe @sprintf("op = %3.1f,  sc1 = %10.4e, sc2 = %10.4e, SF1[end] = %8.2e, SF2[end] = %8.2e, Load/y = %8.2e", op, stg1.SC, stg2.SC, stg1.SF[end], stg2.SF[end], Load_per_year)
        @infoe @sprintf("IF1 = %10.4e, OF1 = %10.4e, IF1-OF1 = %10.4e, SF1[end] = %10.4e, CT1 = %10.4e, IM1 = %10.4e", IF1, OF1,  IF1-OF1, SF1, CT1, IM1)
        @infoe @sprintf("IF2 = %10.4e, OF2 = %10.4e, IF2-OF2 = %10.4e, SF2[end] = %10.4e, CT2 = %10.4e, IM2 = %10.4e", IF2, OF2,  IF1-OF2, SF2, CT2, IM2)
    end

    if plot_p
        fig = [1]
        plot_powers(dates, Load_ec, Load, RP_ec, RP, 0, fig_dir, punit, fig)
        plot_detrended(dates, RP, RP_de, RP_trend, ΔEL, Load, Load_de, Load_trend, punit, fig_dir, fig, data_are_averaged = false)
        plot_storage_fill_level(dates, Load_de, RP_de, res1, oprod, "S1", fig_dir, fig, "storage_fill1", punit, plot_all_p = plot_all_p)
        plot_storage_fill_level(dates, Load_de, RP_de, res2, oprod, "S2", fig_dir, fig, "storage_fill2", punit, plot_all_p = plot_all_p)
    end
end

"""
    load data and compute and plot storage fill levels, data are smoothed using moving averages
"""
function comp_and_plot_averaged(stc1, stc2, oprod, data_dir, fig_dir, punit, start_year, stop_year; plot_p = false, plot_all_p = false, do_log = false)
    dates, Load, RP_ec, RP, RP_de, RP_trend, ΔEL, Load_de, Load_trend = load_data(data_dir, punit, start_year, stop_year)

    averaging_hours = 24*7;
    Load_av, dates_av = averaging(Load, dates, averaging_hours, method = "moving_average")
    RP_av,   dates_av = averaging(RP,   dates, averaging_hours, method = "moving_average")

    RP_av_de, RP_av_trend, ΔEL_av, Load_av_de, Load_av_trend = scale_and_detrend(Load_av, RP_av);
    (res1, res2, oprod) = compute_storage_fill_level(dates, Load_ac_de, RP_av_de, punit, stc1, stc2, oprod, do_log = do_log)

    for (stg1, stg2, op) in zip(res1, res2, oprod)
        Load_per_year = powers_to_energy_per_year(dates, Load)
        IF1  = powers_to_energy_per_year(dates, stg1.IF) / Load_per_year
        OF1  = powers_to_energy_per_year(dates, stg1.OF) / Load_per_year
        CT1  = powers_to_energy_per_year(dates, stg1.CT) / Load_per_year
        IM1  = powers_to_energy_per_year(dates, stg1.IM) / Load_per_year
        IF2  = powers_to_energy_per_year(dates, stg2.IF) / Load_per_year
        OF2  = powers_to_energy_per_year(dates, stg2.OF) / Load_per_year
        CT2  = powers_to_energy_per_year(dates, stg2.CT) / Load_per_year
        IM2  = powers_to_energy_per_year(dates, stg2.IM) / Load_per_year
        SF1  = stg1.SF[end] / Load_per_year
        SF2  = stg2.SF[end] / Load_per_year

        @infoe @sprintf("==== comp_and_plot ================================================================================")
        @infoe @sprintf("op = %3.1f,  sc1 = %10.4e, sc2 = %10.4e, SF1[end] = %8.2e, SF2[end] = %8.2e, Load/y = %8.2e", op, stg1.SC, stg2.SC, stg1.SF[end], stg2.SF[end], Load_per_year)
        @infoe @sprintf("IF1 = %10.4e, OF1 = %10.4e, IF1-OF1 = %10.4e, SF1[end] = %10.4e, CT1 = %10.4e, IM1 = %10.4e", IF1, OF1,  IF1-OF1, SF1, CT1, IM1)
        @infoe @sprintf("IF2 = %10.4e, OF2 = %10.4e, IF2-OF2 = %10.4e, SF2[end] = %10.4e, CT2 = %10.4e, IM2 = %10.4e", IF2, OF2,  IF1-OF2, SF2, CT2, IM2)
    end

    if plot_p
        fig_dir = joinpath(fig_dir, "averaged")
        mkpath(fig_dir)
        fig = [1]
        plot_powers(dates, Load_ec, Load_av, RP_ec, RP_av, averaging_hours, fig_dir, punit, fig)
        plot_detrended(dates, RP, RP_av_de, RP_av_trend, ΔEL_av, Load_av, Load_av_de, Load_av_trend, punit, fig_dir, fig, data_are_averaged = true)
        plot_storage_fill_level(dates, Load_av_de, RP_av_de, res1, oprod, "S1", fig_dir, fig, "storage_fill1", punit, plot_all_p = plot_all_p)
        plot_storage_fill_level(dates, Load_av_de, RP_av_de, res2, oprod, "S2", fig_dir, fig, "storage_fill2", punit, plot_all_p = plot_all_p)
    end
end
