using Statistics

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("storage/utils.jl")
include("storage/hdf5_utils.jl")
include("storage/data_energy.jl")

root = dirname(@__DIR__)
get_data_dir()   = joinpath(root, "data_ise", "energy_charts_data")
get_output_dir() = joinpath(root, "data_ise", "output")
get_fig_dir()    = joinpath(root, "data_ise", "figures")


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
    @infoe @sprintf("(mean(RP) / mean(Load)) = %3.2f", 1.0/scale)

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

"""
    compute_storage_level(dates, Load, RP, eunit, over_production, storage_capacity)

    given Load, RP, over_production and storage_capacity compute storage level as a funtion of time

    dates : times
    Load : power consumed
    RP : renewable power production
    eunit : unit of Load and RP (MW, GW, TW)
    over_production : renewable over production capacity factor, 1.0 is no over production capacity
    storage_capacity : storage capacity
"""
function compute_storage_level(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64},
    eunit::String, over_production::Float64, storage_capacity::Float64)

    M2 = uconvert(eunit, "TW")
    Load = @. Load * M2
    RP = @. RP * M2

    # Dates.value(DateTime) => ms since 1 AD
    ms = Dates.value(dates[2] - dates[1])
    # energy-chart data => 15 min time resolution, 1h = 60 * 60 * 1000 ms = 3.6e6 ms
    Δh = ms/(3.6e6)

    # TWh => TW 15 min
    storage_capacity = storage_capacity / Δh

    # renewable energy over production
    RP = RP .* over_production
    # surplus or shortfall power production
    ΔP = (RP - Load)

    n  = size(RP, 1)
    inflow  = zeros(Float64, n)
    outflow = zeros(Float64, n)
    Import  = zeros(Float64, n)
    Export  = zeros(Float64, n)
    storage_fill = ones(Float64, n)

    storage_fill[1] = max(0.0, ΔP[1])
    for i in 2:n
        # if surplus
        if ΔP[i] > 0.0
            # add surplus energy to storage (time intervall is multiplied at the end)
            storage_fill[i] = storage_fill[i-1] + ΔP[i]
            inflow[i] = ΔP[i]
            # if storage_capacity is reached
            if storage_fill[i] > storage_capacity
                D = storage_fill[i] - storage_capacity
                storage_fill[i] = storage_capacity
                inflow[i] = ΔP[i] - D
                Export[i] = D
            end
        else # if shortfall
            # substract shortfall energy from storage (time intervall is multiplied at the end)
            storage_fill[i] = storage_fill[i-1] + ΔP[i]
            outflow[i] = ΔP[i]
            # if stotage is empty
            if storage_fill[i] < 0.0
                D = storage_fill[i]
                storage_fill[i] = 0.0
                outflow[i] = ΔP[i] - D
                Import[i] = D
            end
        end
    end
    storage_fill = storage_fill .* Δh
    (storage_fill, inflow ./ M2, outflow ./ M2, Import ./ M2, Export ./ M2)
end

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function compute_storage_fill_level(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64}, eunit::String)
    # TWh
    storage_capacity = [14.0, 26.0, 35.0, 45.0]
    over_production  = [1.5, 1.2, 1.15, 1.1, 1.05]
    res = []
    for (s,o) in zip(storage_capacity, over_production)
        storage_fill, inflow, outflow, Import, Export = compute_storage_level(dates, Load, RP, eunit, o, s)
        push!(res, (storage_fill, inflow, outflow, Import, Export, o, s))
    end
    res
end

"""
    determine_overproduction(Load, RP)

    determine minimum storage capacity as a function of over production

    Load : detrended Load
    RP : detrended and scaled renewable power
    eunit : unit of Load and RP (MW. GW, TW)
"""
function determine_overproduction(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64}, eunit::String)
    overproduction = collect(LinRange(1.05, 1.5, 20))
    storage_capacities = []
    for op in overproduction
        storage_capacity = 1.0
        minS = -1.0
        it = 0
        while minS < 0.0 && it < 50
            storage_level = compute_storage_level(dates, Load, RP, eunit, op, storage_capacity)
            min_storage_level = minimum(storage_level)
            storage_capacity = storage_capacity - min_storage_level
            it += 1
        end
        push!(storage_capacities, storage_capacity)
    end
    plt.plot(overproduction, storage_capacities, "r.")
end

"""
    plot_powers(dates, Load, RP, averaging_hours, fig_dir, eunit, fig)

    plot Load and RP

    dates - times
    Load -
    RP   - 
    averaging_hours - if data are averaged number of hours to avergae over
    fig_dir - directory where figures are saved to
    eunit - (MW, GW, TW)
    fig - number of matplotlib figure
"""
function plot_powers(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64},
    averaging_hours::Int64, fig_dir::String, eunit::String, fig::Vector{Int64})

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, Load, label="Load")
    pl.plot(dates, RP, label="Renewables")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", eunit))
    pl.grid()
    pl.title("energy_charts data")
    pl.legend()
    if averaging_hours > 1
        pl.title(@sprintf("averaged, %d days window", div(averaging_hours,24)))
        pl.savefig(joinpath(fig_dir, @sprintf("RP_averaged.png")))
    else
        pl.savefig(joinpath(fig_dir, @sprintf("RP.png")))
    end
end

function plot_detrended(dates::Vector{DateTime}, RP::Vector{Float64}, RP_de::Vector{Float64}, RP_trend::Vector{Float64}, ΔEL::Vector{Float64},
    Load::Vector{Float64}, Load_de::Vector{Float64}, Load_trend::Vector{Float64}, 
    eunit ::String, fig_dir::String, fig::Vector{Int64}; data_are_averaged = false)

    eunit_factor = uconvert(eunit, "GW")

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
    pl.plot(dates, Load        .* eunit_factor, label = label1)
    pl.plot(dates, Load_de     .* eunit_factor, label = label2)
    pl.plot(dates, Load_trend  .* eunit_factor, "r", linewidth=3, label = "Load_trend")
    pl.xlabel("time")
    pl.ylabel("P [GW]")
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
    pl.plot(dates, RP       .* eunit_factor, label = label1)
    pl.plot(dates, RP_de    .* eunit_factor, label = label2)
    pl.plot(dates, RP_trend .* eunit_factor,  "r", linewidth=3, label = "RP_trend")
    pl.xlabel("time")
    pl.ylabel("P [GW]")
    pl.grid()
    pl.legend()
    pl.title("Renewable detrended")
    pl.savefig(joinpath(fig_dir, path1))

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, ΔEL .* eunit_factor, label = label3)
    pl.xlabel("time")
    pl.ylabel("(P_E - P_L) [GW]")
    pl.legend()
    pl.grid()
    pl.title("Renewable power - Load")
    pl.savefig(joinpath(fig_dir, path2))
end

function plot_storage_fill_level(dates::Vector{DateTime}, Load::Vector{Float64}, RP::Vector{Float64}, 
    storage_fill_res, fig_dir::String, fig::Vector{Int64}, pngpath::String; plot_all_p = false)

    prozent = "%"
    pl.figure(fig[1]); fig[1] += 1
    for (storage_fill, inflow, outflow, Import, Export, o, s) in storage_fill_res
        label = @sprintf("storage_cpacity=%2.f TWh, op = %3.f %s", s, (o-1.0)*100.0, prozent)
        pl.plot(dates, storage_fill, label = label)
    end
    pl.xlabel("time")
    pl.ylabel("storage fill level [TWh]")
    pl.legend()
    pl.grid()
    pl.savefig(joinpath(fig_dir, pngpath.*".png"))

    if plot_all_p
        for (storage_fill, inflow, outflow, Import, Export, o, s) in storage_fill_res
            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, RP .*o ,   label="RP")
            pl.plot(dates, Load,   label="Load")
            pl.title(@sprintf("Load, RP, %2.f TWh, %3.f %s", s, (o-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("1_%2.f_%3.f.png", s, (o-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, inflow,   label="inflow")
            pl.plot(dates, Load,   label="Load")
            pl.title(@sprintf("Storage power inflow %2.f TWh, %3.f %s", s, (o-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("2_%2.f_%3.f.png", s, (o-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, outflow,  label="outflow")
            pl.title(@sprintf("Storage power outflow %2.f TWh, %3.f %s", s, (o-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("3_%2.f_%3.f.png", s, (o-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Import, label="import")
            pl.title(@sprintf("Power import %2.f TWh, %3.f %s", s, (o-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("4_%2.f_%3.f.png", s, (o-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, Export, label="export")
            pl.title(@sprintf("Power export %2.f TWh, %3.f %s", s, (o-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("5_%2.f_%3.f.png", s, (o-1.0)*100.0)))
        end
    end
end

"""
    load data and compute and plot storage fille levels, original times (15 min)
"""
function comp_and_plot(;plot_p = false)
    fig_dir = get_fig_dir()

    # use GW
    eunit = ["MW", "GW", "TW"][2]
    data = EnergyData(get_data_dir(), eunit);
    Load = data.Load
    # renewables is sum of wind onshore, wind offshore and solar
    RP = @. data.Woff + data.Won + data.Solar;
    dates = data.dates

    RP_de, RP_trend, ΔEL, Load_de, Load_trend = scale_and_detrend(Load, RP);
    storage_fill_res = compute_storage_fill_level(dates, Load_de, RP_de, eunit)

    if plot_p
        fig = [1]
        plot_powers(dates, Load, RP, 0, fig_dir, eunit, fig)
        plot_detrended(dates, RP, RP_de, RP_trend, ΔEL, Load, Load_de, Load_trend, eunit, fig_dir, fig, data_are_averaged = false)
        plot_storage_fill_level(dates, Load_de, RP_de, storage_fill_res, fig_dir, fig, "storage_fill")
    end
end

"""
    load data and compute and plot storage fill levels, data are smoothed using moving averages
"""
function comp_and_plot_averaged(;plot_p = false)
    fig_dir = get_fig_dir()

    # use GW
    eunit = ["MW", "GW", "TW"][2]
    data = EnergyData(get_data_dir(), eunit);
    Load = data.Load
    RP = @. data.Woff + data.Won + data.Solar;
    dates = data.dates

    averaging_hours = 24*7;
    Load_av, dates_av = averaging(Load, dates, averaging_hours, method = "moving_average");
    RP_av,   dates_av = averaging(RP,   dates, averaging_hours, method = "moving_average");

    RP_av_de, RP_av_trend, ΔEL_av, Load_av_de, Load_av_trend = scale_and_detrend(Load_av, RP_av);
    storage_fill_res = compute_storage_fill_level(dates, Load_av_de, RP_av_de, eunit)

    if plot_p
        fig_dir = joinpath(fig_dir, "averaged")
        mkpath(fig_dir)
        fig = [1]
        plot_powers(dates, Load_av, RP_av, averaging_hours, fig_dir, eunit, fig)
        plot_detrended(dates, RP, RP_av_de, RP_av_trend, ΔEL_av, Load_av, Load_av_de, Load_av_trend, eunit, fig_dir, fig, data_are_averaged = true)
        plot_storage_fill_level(dates, Load_av_de, RP_av_de, storage_fill_res, fig_dir, fig, "storage_file_av")
    end
end

comp_and_plot(plot_p=true);

comp_and_plot_averaged(plot_p=true);