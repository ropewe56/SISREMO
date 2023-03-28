using Common
using Statistics

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
    dates
    returns ΔTh - multiple of an hour per step
"""
function get_step_ΔTh(dates)
    # Dates.value(DateTime) => ms since 1 AD
    ms = Dates.value(dates[2] - dates[1])
    # energy-chart data => 15 min time resolution ΔTh = 0.25
    # energy-chart data averaged over 4 steps => 1h time resolution ΔTh = 1.0
    # 1h = 60 * 60 * 1000 ms = 3.6e6 ms
    ms/(3.6e6)
end

function nb_years(dates)
    ΔT_ms_e  = Dates.value(dates[end] - dates[1])
    Δh_e     = ΔT_ms_e/(3.6e6)
    Δh_e/(365*24)
end
function powers_to_energy_per_year(dates, P)
    ΔTh = get_step_ΔTh(dates)
    nb_years = nb_years(dates)
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
    SC :: Float64         # storage capacity [punit_h]
    SF :: Vector{Float64} # storage fill     [punit_h]
    IF :: Vector{Float64} # inflow           [punit]
    OF :: Vector{Float64} # outflow          [punit]
    CT :: Vector{Float64} # curtailement     [punit]
    IM :: Vector{Float64} # import           [punit]
end
function Storage(SC::Float64, n::Int64)
    SF = zeros(Float64, n)
    IF = zeros(Float64, n)
    OF = zeros(Float64, n)
    CT = zeros(Float64, n)
    IM = zeros(Float64, n)
    Storage(SC, SF, IF, OF, CT, IM)
end

function write_to_log(stg::Storage, ΔPin, out, i, j)
    B = stg.IF[i] - stg.OF[i] - (stg.SF[i] - stg.SF[i-1])
    write(out, @sprintf("%5d  %d  ΔP = %9.2e, S = %9.2e, I = %9.2e, O = %9.2e, C = %9.2e, M = %9.2e, B = %9.2e\n", i, j, ΔPin, stg.SF[i], stg.IF[i], stg.OF[i], stg.CT[i], stg.IM[i], B))
end

function charge(stg::Storage, ΔPin, i; out = nothing)
    local ΔPout
    ΔSC = stg.SC - stg.SF[i-1]
    if ΔPin > 0.0
        if ΔSC > ΔPin
            stg.IF[i] = ΔPin
            stg.SF[i] = stg.SF[i-1] + stg.IF[i]
            ΔPout    = 0.0     # Load is satisfied

            if out !== nothing write_to_log(stg, ΔPin, out, i, 1) end
        else
            stg.IF[i] = ΔSC
            stg.SF[i] = stg.SF[i-1] + stg.IF[i]
            stg.CT[i] = ΔPin - stg.IF[i] # to be curtailed (or exported)
            ΔPout     = stg.CT[i] # Load is staisfied, ΔPout remains

            if out !== nothing write_to_log(stg, ΔPin, out, i, 2) end
        end
    else
        if stg.SF[i-1] > -ΔPin   # ΔPin < 0
            stg.OF[i] = -ΔPin    # stg.OF[i] < 0
            stg.SF[i] = stg.SF[i-1] - stg.OF[i]
            ΔPout     = 0.0       # Load is satisfied

            if out !== nothing write_to_log(stg, ΔPin, out, i, 3) end
        else
            stg.OF[i] = stg.SF[i-1]
            stg.SF[i] = stg.SF[i-1] - stg.OF[i]
            stg.IM[i] = -ΔPin - stg.OF[i] # to be imported
            ΔPout     = -stg.IM[i]        # missing for Load to be satisfied

            if out !== nothing write_to_log(stg, ΔPin, out, i, 4) end
        end
    end
    ΔPout
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
    stc1::Float64, stc2::Float64, oprod::Float64; do_log = false)

    # renewable energy over production
    RP = RP .* oprod
    # surplus or shortfall power production
    ΔP = (RP - Load)

    out1, out2, out3 = nothing, nothing, nothing
    if do_log
        out1 = open(@sprintf("log_%3.1f_1.log", oprod), "w")
        out2 = open(@sprintf("log_%3.1f_2.log", oprod), "w")
        out3 = open(@sprintf("log_%3.1f_3.log", oprod), "w")
    end

    n  = size(RP, 1)
    stg1 = Storage(stc1, n)
    stg2 = Storage(stc2, n)
    for i in 2:n
        ΔPout1 = charge(stg1, ΔP[i],  i, out = out1)
        ΔPout2 = charge(stg2, ΔPout1, i, out = out2)

        balance1 = stg1.IF[i] - stg1.OF[i] - (stg1.SF[i] - stg1.SF[i-1])
        balance2 = stg2.IF[i] - stg2.OF[i] - (stg2.SF[i] - stg2.SF[i-1])
        if abs(balance1) > 1.0e-10 || abs(balance2) > 1.0e-10
            write(out3, @sprintf("%d  b1 = %9.2e  b2 = %9.2e  P1 = %9.2e  P2 = %9.2e\n", i, balance1, balance2, ΔPout1, ΔPout2))
        end
    end

    @infoe @sprintf("op = %3.1f, SC1 = %8.2e, SC2 = %8.2e, SF1[end] = %8.2e, SF2[end] = %8.2e", oprod, stg1.SC, stg2.SC, stg1.SF[end], stg2.SF[end])

    IF1 = sum(stg1.IF)
    OF1 = sum(stg1.OF)
    SF1 = stg1.SF[end]

    IF2 = sum(stg2.IF)
    OF2 = sum(stg2.OF)
    SF2 = stg2.SF[end]

    @infoe @sprintf("IF1 = %9.3e, OF1 = %9.3e, SF1[end] = %9.3e, B1 = %9.3e", IF1, OF1, SF1,  IF1 - OF1 - SF1)
    @infoe @sprintf("IF2 = %9.3e, OF2 = %9.3e, SF2[end] = %9.3e, B2 = %9.3e", IF2, OF2, SF2,  IF2 - OF2 - SF2)

    (stg1, stg2)
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
function plot_powers(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64},
    averaging_hours::Int64, fig_dir::String, punit::String, fig::Vector{Int64})

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, Load, label="Load")
    pl.plot(dates, RP, label="Renewables")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", punit))
    pl.grid()
    pl.title("energy_charts data")
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
    pl.ylabel(@sprintf("storage fill level [%s]", punit))
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
            pl.savefig(joinpath(fig_dir, @sprintf("%s_1_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Load,   label="Load")
            pl.plot(dates, stg.IF,   label="IF")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Storage power inflow: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_2_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Load,   label="Load")
            pl.plot(dates, stg.OF,  label="OF")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Storage power outflow: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_3_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Load,   label="Load")
            pl.plot(dates, stg.IM, label="IM")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Power import: SC=%2.f[%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_4_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Load,   label="Load")
            pl.plot(dates, stg.CT, label="CT")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Power curtailment: SC=%2.f %s, op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_5_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))
        end
    end
end

function load_data(data_dir, punit)
    data = EnergyData(data_dir, punit);
    Load = data.Load

    # renewables is sum of wind onshore, wind offshore and solar
    RP = @. data.Woff + data.Won + data.Solar;

    # Biomass
    BM = data.Bio

    # dates
    dates = data.dates
    @infoe @sprintf("start: (%d, %d, %d), stop: (%d, %d, %d)",
        Dates.day(dates[1]),   Dates.month(dates[1]),   Dates.year(dates[1]),
        Dates.day(dates[end]), Dates.month(dates[end]), Dates.year(dates[end]))

    # energy-chart data => 15 min time resolution, 1h = 60 * 60 * 1000 ms = 3.6e6 ms
    @infoe @sprintf("RE = %8.2e %sh", energy(dates, RP  ), punit)
    @infoe @sprintf("LE = %8.2e %sh", energy(dates, Load), punit)

    RP_de, RP_trend, ΔEL, Load_de, Load_trend = scale_and_detrend(Load, RP);

    @infoe @sprintf("RE_de = %8.2e %sh", energy(dates, RP_de  ), punit)
    @infoe @sprintf("LE_de = %8.2e %sh", energy(dates, Load_de), punit)

    dates, Load, RP, RP_de, RP_trend, ΔEL, Load_de, Load_trend, BM
end

"""
    load data and compute and plot storage fille levels, original times (15 min)
"""
function comp_and_plot(stc1, stc2, oprod, data_dir, fig_dir, punit; plot_p = false, plot_all_p = false, do_log = false)

    dates, Load, RP, RP_de, RP_trend, ΔEL, Load_de, Load_trend, BM = load_data(data_dir, punit)
    RP_de = RP_de + BM

    (res1, res2, oprod) = compute_storage_fill_level(dates, Load_de, RP_de, punit, stc1, stc2, oprod, do_log = do_log)
    # Vector{(storage_fill, stg1, stg2, sc, op)}

    for (stg1, stg2, op) in zip(res1, res2, oprod)
        Load_per_year = energy(dates, Load)
        IF1  = energy(dates, stg1.IF) / Load_per_year
        OF1  = energy(dates, stg1.OF) / Load_per_year
        CT1  = energy(dates, stg1.CT) / Load_per_year
        IM1  = energy(dates, stg1.IM) / Load_per_year
        IF2  = energy(dates, stg2.IF) / Load_per_year
        OF2  = energy(dates, stg2.OF) / Load_per_year
        CT2  = energy(dates, stg2.CT) / Load_per_year
        IM2  = energy(dates, stg2.IM) / Load_per_year
        SF1  = stg1.SF[end] / Load_per_year
        SF2  = stg2.SF[end] / Load_per_year

        @infoe @sprintf("==== comp_and_plot ================================================================================")
        @infoe @sprintf("op = %3.1f,  sc1 = %10.4e, sc2 = %10.4e, SF1[end] = %8.2e, SF2[end] = %8.2e, Load/y = %8.2e", op, stg1.SC, stg2.SC, stg1.SF[end], stg2.SF[end], Load_per_year)
        @infoe @sprintf("IF1 = %10.4e, OF1 = %10.4e, IF1-OF1 = %10.4e, SF1[end] = %10.4e, CT1 = %10.4e, IM1 = %10.4e", IF1, OF1,  IF1-OF1, SF1, CT1, IM1)
        @infoe @sprintf("IF2 = %10.4e, OF2 = %10.4e, IF2-OF2 = %10.4e, SF2[end] = %10.4e, CT2 = %10.4e, IM2 = %10.4e", IF2, OF2,  IF1-OF2, SF2, CT2, IM2)

    end

    if plot_p
        fig = [1]
        plot_powers(dates, Load, RP, 0, fig_dir, punit, fig)
        plot_detrended(dates, RP, RP_de, RP_trend, ΔEL, Load, Load_de, Load_trend, punit, fig_dir, fig, data_are_averaged = false)
        plot_storage_fill_level(dates, Load_de, RP_de, res1, oprod, "S1", fig_dir, fig, "storage_fill1", eunit, plot_all_p = plot_all_p)
        plot_storage_fill_level(dates, Load_de, RP_de, res2, oprod, "S2", fig_dir, fig, "storage_fill2", eunit, plot_all_p = plot_all_p)
    end
end

"""
    load data and compute and plot storage fill levels, data are smoothed using moving averages
"""
function comp_and_plot_averaged(data_dir, fig_dir, punit; plot_p = false)

    data = EnergyData(data_dir, punit);
    Load = data.Load
    RP = @. data.Woff + data.Won + data.Solar;
    dates = data.dates

    averaging_hours = 24*7;
    Load_av, dates_av = averaging(Load, dates, averaging_hours, method = "moving_average");
    RP_av,   dates_av = averaging(RP,   dates, averaging_hours, method = "moving_average");

    RP_av_de, RP_av_trend, ΔEL_av, Load_av_de, Load_av_trend = scale_and_detrend(Load_av, RP_av);
    storage_fill_res = compute_storage_fill_level(dates, Load_av_de, RP_av_de, punit)

    if plot_p
        fig_dir = joinpath(fig_dir, "averaged")
        mkpath(fig_dir)
        fig = [1]
        plot_powers(dates, Load_av, RP_av, averaging_hours, fig_dir, punit, fig)
        plot_detrended(dates, RP, RP_av_de, RP_av_trend, ΔEL_av, Load_av, Load_av_de, Load_av_trend, punit, fig_dir, fig, data_are_averaged = true)
        plot_storage_fill_level(dates, Load_av_de, RP_av_de, storage_fill_res, fig_dir, fig, "storage_file_av", punit)
    end
end
