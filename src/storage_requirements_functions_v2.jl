using Common
using Statistics
using Interpolations

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("storage/utils.jl")
#include("storage/hdf5_utils.jl")
include("storage/data_energy.jl")

include("plot_results_v2.jl")

root = dirname(@__DIR__)
get_data_dir()   = joinpath(root, "data")
get_fig_dir()    = joinpath(root, "figures")
mkpath(get_fig_dir())

function detrend_times_series(A)
    n  = length(A)
    nh = div(n,2)
    k = 2
    A_trend = polynomial_fit(A, k)
    A_de = @. A * A_trend[nh] / A_trend
    A_de, A_trend
end

function renewables_scale_and_detrend(data_dir, data, punit, start_year, stop_year)
    IP = InstalledPower(data, data_dir, punit, start_year, stop_year)

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

    s1 = mean(data.Woff) / mean(Woff)
    s2 = mean(data.Won)  / mean(Won)
    s3 = mean(data.Solar)/ mean(Solar)
    s4 = mean(data.Bio)  / mean(Bio)

    Woff  = detrend_time_series(Woff  .* s1)
    Won   = detrend_time_series(Won   .* s2)
    Solar = detrend_time_series(Solar .* s3)
    Bio   = detrend_time_series(Bio   .* s4)

    L_de = detrend_time_series(data.Load)
    Dict("Load" => L_de, "Woff" => Woff .* s1, "Wo" => Won .* s2, "Solar" => Solar .* s3, "Bio" => Bio .* s4)
end

function renewables_detrend(data)
    Woff  = detrend_time_series(data.Woff )
    Won   = detrend_time_series(data.Won  )
    Solar = detrend_time_series(data.Solar)
    Bio   = detrend_time_series(data.Bio  )

    L_de = detrend_time_series(data.Load)
    Dict("Load" => L_de, "Woff" => Woff, "Wo" => Won, "Solar" => Solar, "Bio" => Bio)
end


function scale_detrend_print_info(data_ec, data)

    # dates
    dates = data_ec.dates
    @infoe @sprintf("start: (%d, %d, %d), stop: (%d, %d, %d)",
        Dates.day(dates[1]),   Dates.month(dates[1]),   Dates.year(dates[1]),
        Dates.day(dates[end]), Dates.month(dates[end]), Dates.year(dates[end]))

    # energy-chart data => 15 min time resolution, 1h = 60 * 60 * 1000 ms = 3.6e6 ms
    @infoe @sprintf("PE = %8.2e %sh", powers_to_energy_per_year(dates, P_ec), punit)
    @infoe @sprintf("LE = %8.2e %sh", powers_to_energy_per_year(dates, L_ec), punit)


    @infoe @sprintf("PE_de = %8.2e %sh", powers_to_energy_per_year(dates, P_de  ), punit)
    @infoe @sprintf("LE_de = %8.2e %sh", powers_to_energy_per_year(dates, L_de), punit)

    dates, L_ec, L, P_ec, P, P_de, P_trend, ΔEL, L_de, L_trend
end

function load_data(data_dir, punit, start_year, stop_year; scale_to_installed_power = true, plot_p=false)

    data = PowerData(data_dir, punit, start_year, stop_year);
    L_ec = data.Load
    dates = data.dates

    # energy charts
    P_ec = @. data.Woff + data.Won + data.Solar + data.Bio;


    # dates
    @infoe @sprintf("start: (%d, %d, %d), stop: (%d, %d, %d)",
        Dates.day(dates[1]),   Dates.month(dates[1]),   Dates.year(dates[1]),
        Dates.day(dates[end]), Dates.month(dates[end]), Dates.year(dates[end]))

    # energy-chart data => 15 min time resolution, 1h = 60 * 60 * 1000 ms = 3.6e6 ms
    @infoe @sprintf("PE = %8.2e %sh", powers_to_energy_per_year(dates, P_sc), punit)
    @infoe @sprintf("LE = %8.2e %sh", powers_to_energy_per_year(dates, L_sc), punit)

    P_de, P_trend, ΔEL, L_de, L_trend = scale_and_detrend(L_ec, P_ec);

    @infoe @sprintf("PE_de = %8.2e %sh", powers_to_energy_per_year(dates, P_de), punit)
    @infoe @sprintf("LE_de = %8.2e %sh", powers_to_energy_per_year(dates, L_de), punit)


    dates, data.Load, L, P_ec, P, P_de, P_trend, ΔEL, L_de, L_trend
end

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

function power_step(stg::Storage, L, P, i)
    ΔSC = stg.SC - stg.SF[i-1]
    k = 0
    if P - L > 0.0
        stg.I2[i] = L
        stg.I3[i] = P - L
        stg.I4[i] = min(stg.I3[i], ΔSC/stg.ηin)
        stg.SF[i] = stg.SF[i-1] + stg.I4[i]*stg.ηin
        stg.I5[i] = max(0.0, stg.I3[i]-stg.I4[i])  # curtailment
        k = 1
    else
        stg.I2[i] = P
        stg.I6[i] = min(stg.SF[i-1] * stg.ηout, L - P)
        stg.SF[i] = stg.SF[i-1] - stg.I6[i]/stg.ηout
        stg.I7[i] = L - stg.I2[i] - stg.I6[i]      # residual load
        k = 2
    end
    k
end

function power_curtailment(powers, op)
    L  = powers["Load"]
    B  = powers["Bio"]
    Wf = powers["Woff"]
    Wn = powers["Won"]
    So = powers["Solar"]
    R  = Wf .+ Wn .+ So

    mean_L = mean(L)
    mean_B = mean(B)
    mean_R = mean(R)

    # mean_L*op = mean_R * scale + mean_B
    scale = (mean_L*op - mean_B) / mean_R
    P = R .* scale .+ Bio

    L, P
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
function compute_storage_level(powers, st_capacities, op::Float64; log_p=false)

    # curtailment
    L, P = power_curtailment(powers, op)

    storages = []
    nb_steps  = size(P, 1)
    for stc in st_capacities
        storage = Storage(stc, nb_steps)
        storage.SF[1] = stc*0.5
        push!(storages, storage)
    end

    if log_p
        out1 = open(joinpath(@__DIR__, "log1.log"), "w")
        out2 = open(joinpath(@__DIR__, "log2.log"), "w")
        for i in 2:nb_steps
            k1 = power_step(storages[1], L[i], P[i], i)

            write(out1, @sprintf("%3d, %d, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n", i, k1, L[i], P[i],
                storages[1].I2[i],
                storages[1].I3[i],
                storages[1].I4[i],
                storages[1].I5[i],
                storages[1].I6[i],
                storages[1].I7[i],
                storages[1].SF[i],
                ))
            for j in 2:length(storages)
                k2 = power_step(storages[j], storages[j-1].I7[i], storages[j-1].I5[i], i)

                write(out2, @sprintf("%3d, %d, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n", i, k2, storages[j-1].I7[i], storages[j-1].I5[i],
                storages[j].I2[i],
                storages[j].I3[i],
                storages[j].I4[i],
                storages[j].I5[i],
                storages[j].I6[i],
                storages[j].I7[i],
                storages[j].SF[i],
                ))
            end
        end
        close(out1)
        close(out2)
    else
        for i in 2:nb_steps
            k1 = power_step(storages[1], L[i], P[i], i)
            for j in 2:length(storages)
                k2 = power_step(storages[j], storages[j-1].I7[i], storages[j-1].I5[i], i)
            end
        end
    end

    sum_I2 = 0.0
    sum_I6 = 0.0
    for (j, stg) in enumerate(storages)
        direct = sum(stg.I2)
        st_in  = sum(stg.I4)
        st_out = sum(stg.I6)
        other  = sum(stg.I7)
        curt   = sum(stg.I5)
        @infoe @sprintf("j = %d, op = %3.1f, SC = %8.2e, SF[end] = %8.2e", j, oprod, stg.SC, stg.SF[end])
        @infoe @sprintf("direct = %8.2e, st_in = %8.2e, st_out = %8.2e, other = %8.2e, curt = %8.2e", direct, st_in, st_out, other, curt)

        sum_I2 += direct
        sum_I6 += st_out
    end

    sum_I7 = sum(storages[end].I7)
    sum_I5 = sum(storages[end].I5)
    @infoe (sum_I2 + sum_I6 + sum_I7), sum_I5

    storages
end

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function compute_storage_fill_level(powers, st_capacities::Vector{Vector{Float64}}, oprod::Vector{Float64})
    @infoe @sprintf("==== compute_storage_fill_level ===================================================================")

    results = []
    for (i, op) in enumerate(oprod)
        storages = compute_storage_level(powers, st_capacities[i], op)
        push!(results, storages)
    end
    results
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
    load data and compute and plot storage fille levels, original times (15 min)
"""
function comp_and_plot(st_capacities, oprod, data_dir, fig_dir, punit, start_year, stop_year; plot_p = false, plot_all_p = false, do_log = false, scaled_to_installed_power_p = true)

    data_ec = PowerData(data_dir, punit, start_year, stop_year);
    local powers
    if scaled_to_installed_power_p
        powers = renewables_scale_and_detrend(data_dir, data_ec, punit, start_year, stop_year);
    else
        powers = renewables_detrend(data_ec);
    end

    results = compute_storage_fill_level(powers, st_capacities, oprod)

    nb_stg = length(results[1])
    stores = []
    for j in 1:nb_stg
        push!(stores, [])
    end
    for j in 1:nb_stg
        for storages in results
            push!(stores[j], storages[j])
        end
    end

    if plot_p
        fig = [1]
        plot_powers(dates, L_ec, L, P_ec, P, 0, fig_dir, punit, fig)
        plot_detrended(dates, P, P_de, P_trend, ΔEL, L, L_de, L_trend, punit, fig_dir, fig, data_are_averaged = false)

        for j in 1:nb_stg
            plot_storage_fill_level(dates, L_de, P_de, stores[j], oprod, j, fig_dir, fig, punit, plot_all_p = plot_all_p)
        end
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
