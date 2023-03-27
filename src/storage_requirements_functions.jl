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
    CP :: Float64         # storage capacity [punit_h]
    SF :: Vector{Float64} # storage fill     [punit_h]
    IF :: Vector{Float64} # inflow           [punit]
    OF :: Vector{Float64} # outflow          [punit]
    CT :: Vector{Float64} # curtailement     [punit]
    IM :: Vector{Float64} # import           [punit]
end
function Storage(CP::Float64, n::Int64)
    SF = zeros(Float64, n)
    IF = zeros(Float64, n)
    OF = zeros(Float64, n)
    CT = zeros(Float64, n)
    IM = zeros(Float64, n)
    Storage(CP, SF, IF, OF, CT, IM)
end

function charge(st::Storage, ΔPin, i, out, j)
    ΔPout = 0.0
    if ΔPin > 0.0
        D = st.CP - st.SF[i-1]
        if D > ΔPin
            st.SF[i] = st.SF[i-1] + ΔPin
            st.IF[i] = ΔPin
            write(out, @sprintf("%d 1 %9.2e,  %9.2e, %9.2e\n", j, st.SF[i], ΔPin, st.CP))
        else
            st.SF[i] = st.SF[i-1] + D
            st.IF[i] = D
            st.CT[i] = ΔPin - D
            ΔPout = st.CT[i]
            write(out, @sprintf("%d 2 %9.2e,  %9.2e, %9.2e\n", j, st.SF[i], ΔPin, st.CP))
        end
    else
        D = st.SF[i-1]
        if D > -ΔPin
            st.SF[i] = st.SF[i-1] + ΔPin
            st.OF[i] = ΔPin
            write(out, @sprintf("%d 3 %9.2e,  %9.2e, %9.2e\n", j, st.SF[i], ΔPin, st.CP))
        else
            st.SF[i] = st.SF[i-1] + D
            st.OF[i] = D
            st.IM[i] = ΔPin - D
            ΔPout = st.IM[i]
            write(out, @sprintf("%d 4 %9.2e,  %9.2e, %9.2e\n", j, st.SF[i], ΔPin, st.CP))
        end
    end
    ΔPout
end


"""
    compute_storage_level(dates, Load, RP, punit, over_production, storage_capacity)

    given Load, RP, over_production and storage_capacity compute storage level as a funtion of time

    dates : times
    Load  : power consumed
    RP    : renewable power production
    punit : power unit of Load and RP (MW, GW, TW)
    over_production  : renewable over production capacity factor, 1.0 is no over production capacity
    storage_capacity : storage capacity
"""
function compute_storage_level(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64},
    stcap1::Float64, stcap2::Float64, oprod::Float64)

    # Dates.value(DateTime) => ms since 1 AD
    ms = Dates.value(dates[2] - dates[1])
    # energy-chart data => 15 min time resolution, 1h = 60 * 60 * 1000 ms = 3.6e6 ms
    Δh = ms/(3.6e6)

    # E => P * 15 min
    stcap1 = stcap1 / Δh
    stcap2 = stcap2 / Δh

    # renewable energy over production
    RP = RP .* oprod
    # surplus or shortfall power production
    ΔP = (RP - Load)

    n  = size(RP, 1)
    out1 = open("llog1.log", "w")
    out2 = open("llog2.log", "w")

    stg1 = Storage(stcap1, n)
    stg2 = Storage(stcap2, n)
    for i in 2:n
        ΔPout1 = charge(stg1, ΔP[i],  i, out1, 1)
        ΔPout2 = charge(stg2, ΔPout1, i, out2, 2)
    end

    stg1.SF = stg1.SF .* Δh
    stg2.SF = stg2.SF .* Δh

    (stg1, stg2)
end

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function compute_storage_fill_level(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64}, punit::String)
    if punit == "GW"
        c = 1.0e3
    elseif punit == "TW"
        c = 1.0
    end
    stcap1  = [14.0, 26.0, 35.0, 45.0] .* c
    oprod  = [1.5, 1.2, 1.15, 1.1, 1.05]
    stcap2 = @. stcap1 * 0.1
    stcap1 = @. stcap1 - stcap2
    @infoe @sprintf("stcap1 = %s, stcap2 = %s", stcap1, stcap2)

    res1 = []
    res2 = []
    scop = []
    for (sc1, sc2, op) in zip(stcap1, stcap2, oprod)
        stg1, stg2 = compute_storage_level(dates, Load, RP, sc1, sc2, op)

        @infoe @sprintf("sc1 = %8.2e, sc2 = %8.2e, op = %3.1f, SF1[end] = %8.2e, CP1 = %8.2e, SF2[end] = %8.2e, CP2 = %8.2e", 
           sc1, sc2, op, stg1.SF[end], stg1.CP, stg2.SF[end], stg2.CP)

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
        stcap1 = 1.0
        stcap2 = 1.0
        minS = -1.0
        it = 0
        while minS < 0.0 && it < 50
            stg1, stg2 = compute_storage_level(dates, Load, RP, op, stcap1, stcap2)

            min_storage_level1 = minimum(stg1.SF)
            min_storage_level2 = minimum(stg2.SF)

            stcap1 = stcap1 - min_storage_level1
            stcap2 = stcap2 - min_storage_level2
            it += 1
        end
        push!(storage_capacities, (stcap1, stcap2))
    end
    plt.plot(overproduction, storage_capacities, "r.")
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
    storage, sc_op, title, fig_dir::String, fig::Vector{Int64}, pngpath::String, punit ::String; plot_all_p = false)

    prozent = "%"
    pl.figure(fig[1]); fig[1] += 1
    for (stg, so) in zip(storage, sc_op)
        sc = so[1]
        op = so[2]
    label = @sprintf("%s storage_cpacity=%2.f %s, op = %3.f %s", title, sc, punit, (op-1.0)*100.0, prozent)
        pl.plot(dates, stg.SF, label = label)
    end
    pl.xlabel("time")
    pl.ylabel(@sprintf("storage fill level [%s]", punit))
    pl.legend()
    pl.grid()
    pl.savefig(joinpath(fig_dir, pngpath.*".png"))

    if plot_all_p
        for (stg, so) in zip(storage, sc_op)
            sc= so[1]
            op= so[2]

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, RP .* op ,   label="RP")
            pl.plot(dates, Load,   label="Load")
            pl.title(@sprintf("%s Load, RP, %2.f %s, %3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_1_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, stg.IF,   label="IF")
            pl.plot(dates, Load,   label="Load")
            pl.title(@sprintf("%s Storage power inflow %2.f %s, %3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_2_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, stg.OF,  label="OF")
            pl.title(@sprintf("%s Storage power outflow %2.f %s, %3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_3_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, stg.IM, label="IM")
            pl.title(@sprintf("%s Power import %2.f %s, %3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_4_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, stg.CT, label="CT")
            pl.title(@sprintf("%s Power curtailment %2.f %s, %3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_5_%2.f_%3.f.png", title, sc, (op-1.0)*100.0)))
        end
    end
end

function time_step(dates)
    ΔT_ms  = Dates.value(dates[2] - dates[1])
    ΔT_ms/(3.6e6)
end
function nb_years(dates)
    ΔT_ms_e  = Dates.value(dates[end] - dates[1])
    Δh_e     = ΔT_ms_e/(3.6e6)
    Δh_e/(365*24)
end
function energy(dates, P)
    ΔT = time_step(dates)
    nb_y = nb_years(dates)
    E = sum(P)*ΔT/nb_y
    #@infoe ("ΔT =", ΔT, "nb_y =", nb_y, "E =", E)
    E
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
function comp_and_plot(data_dir, fig_dir, punit; plot_p = false, plot_all_p = false)

    dates, Load, RP, RP_de, RP_trend, ΔEL, Load_de, Load_trend, BM = load_data(data_dir, punit)
    RP_de = RP_de + BM

    #@infoe (sum(Load)/8.0,  energy(dates, Load))

    (res1, res2, oprod) = compute_storage_fill_level(dates, Load_de, RP_de, punit)
    # Vector{(storage_fill, stg1, stg2, sc, op)}

    for (stg1, stg2, op) in zip(res1, res2, oprod)
        dEL = 1.0/energy(dates, Load)
        IF1  = energy(dates, stg1.IF)*dEL
        OF1  = energy(dates, stg1.OF)*dEL
        CT1  = energy(dates, stg1.CT)*dEL
        IM1  = energy(dates, stg1.IM)*dEL
        IF2  = energy(dates, stg2.IF)*dEL
        OF2  = energy(dates, stg2.OF)*dEL
        CT2  = energy(dates, stg2.CT)*dEL
        IM2  = energy(dates, stg2.IM)*dEL
        SF1  = stg1.SF[end]*dEL
        SF2  = stg2.SF[end]*dEL

        @infoe @sprintf("----------------------")
        @infoe @sprintf("sc = %10.4e, op = %10.4e", stg1.CP, stg2.CP)
        @infoe @sprintf("SF1[end] = %8.2e, SF2[end] = %8.2e, Load/y = %8.2e", stg1.SF[end], stg2.SF[end], energy(dates, Load))
        @infoe @sprintf("IF = %10.4e, OF = %10.4e, IF-OF = %10.4e, SF = %10.4e,  CT = %10.4e, IM = %10.4e", 
            IF1, -OF1,  IF1+OF1, SF1, CT1, IM1)
        @infoe @sprintf("IF = %10.4e, OF = %10.4e, IF-OF = %10.4e, SF = %10.4e,  CT = %10.4e, IM = %10.4e", 
            IF2, -OF2,  IF1+OF2, SF2, CT2, IM2)

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
