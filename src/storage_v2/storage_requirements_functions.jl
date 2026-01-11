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
function number_of_years(dates)
    ΔT_ms_e  = Dates.value(dates[end] - dates[1])
    Δh_e     = ΔT_ms_e/(3.6e6)
    Δh_e/(365*24)
end

"""
    energy produced per year by power series P
"""
function powers_to_energy_per_year(dates, P)
    ΔTh = get_step_ΔTh(dates)
    nb_years = number_of_years(dates)
    sum(P)*ΔTh/nb_years
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

"""
    power_step(stg::Storage, L, P, i)

    L - Load
    P - power
"""
function power_step(stg::Storage, L, P, istep)
    ΔSC = stg.SC - stg.SF[istep-1]
    k = 0
    if P - L > 0.0
        stg.I2[istep] = L
        stg.I3[istep] = P - L
        stg.I4[istep] = min(stg.I3[istep], ΔSC/stg.ηin)
        stg.SF[istep] = stg.SF[istep-1] + stg.I4[istep]*stg.ηin
        stg.I5[istep] = max(0.0, stg.I3[istep]-stg.I4[istep])  # curtailment
        k = 1
    else
        stg.I2[istep] = P
        stg.I6[istep] = min(stg.SF[istep-1] * stg.ηout, L - P)
        stg.SF[istep] = stg.SF[istep-1] - stg.I6[istep]/stg.ηout
        stg.I7[istep] = L - stg.I2[istep] - stg.I6[istep]      # residual load
        k = 2
    end
    k
end

"""
    get_scaled_renewables(powers, op)

    powers : Load, Woff, Won, Solar, Bio
    op : over production factor

    R = Woff+Won+Solar
    s = (<Load>*op - <Bio>) / <Woff+Won+Solar>

    return (Woff+Won+Solar)/<Woff+Won+Solar> * (<Load>*op - <Bio>) + Bio
"""
function get_scaled_renewables(powers, op)
    renewables = powers.Woff .+ powers.Won .+ powers.Solar
    mean_Load = mean(powers.Load)
    mean_Bio = mean(powers.Bio)
    mean_renewables = mean(renewables)

    # mean_L*op = mean_R * scale + mean_B
    scale = (mean_Load*op - mean_Bio) / mean_renewables

    renewables .* scale .+ powers.Bio
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
function compute_storage_level(load, scaled_renewables, storage_capacities, op::Float64, SF1_factor; log_p=false)

    storages = Vector{Storage}(undef, 0)
    nb_steps  = size(scaled_renewables, 1)
    for stc in storage_capacities
        storage = Storage(stc, nb_steps)
        storage.SF[1] = stc*SF1_factor
        push!(storages, storage)
    end

    if log_p
        out1 = open(joinpath(@__DIR__, "log1.log"), "w")
        out2 = open(joinpath(@__DIR__, "log2.log"), "w")
        for istep in 2:nb_steps
            k1 = power_step(storages[istep], load[istep], scaled_renewables[istep], istep)

            write(out1, @sprintf("%3d, %d, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n", istep, k1, load[v], scaled_renewables[istep],
                storages[1].I2[i],
                storages[1].I3[i],
                storages[1].I4[i],
                storages[1].I5[i],
                storages[1].I6[i],
                storages[1].I7[i],
                storages[1].SF[i],
                ))

            nb_strgs = length(storages)

            for j in 2:nb_strgs

                k2 = power_step(storages[j], storages[j-1].I7[istep], storages[j-1].I5[istep], istep)

                write(out2, @sprintf("%3d, %d, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n", 
                    istep, k2, storages[j-1].I7[istep], storages[j-1].I5[istep],
                storages[j].I2[istep],
                storages[j].I3[istep],
                storages[j].I4[istep],
                storages[j].I5[istep],
                storages[j].I6[istep],
                storages[j].I7[istep],
                storages[j].SF[istep],
                ))
            end
        end
        close(out1)
        close(out2)
    else
        for istep in 2:nb_steps
            k1 = power_step(storages[1], load[istep], scaled_renewables[istep], istep)
            nb_strgs = length(storages)
            for j in 2:nb_strgs
                k2 = power_step(storages[j], storages[j-1].I7[istep], storages[j-1].I5[istep], istep)
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
        @info @sprintf("storage_i = %d, op = %3.1f, capacity = %8.2e, storage_fill[end] = %8.2e", j, op, stg.SC, stg.SF[end])
        @info @sprintf("direct = %8.2e, st_in = %8.2e, st_out = %8.2e, other = %8.2e, curt = %8.2e", direct, st_in, st_out, other, curt)

        sum_I2 += direct
        sum_I6 += st_out
    end

    sum_I7 = sum(storages[end].I7)
    sum_I5 = sum(storages[end].I5)
    @info @sprintf("sum: (direct+out+other) = %10.4e, curtailment = %10.4e", (sum_I2 + sum_I6 + sum_I7), sum_I5)

    storages, scaled_renewables
end

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function compute_storage_fill_level(powers, storage_capacities::Vector{Float64}, over_production::Vector{Float64}, SF1_factor)
    @info @sprintf("==== compute_storage_fill_level ===================================================================")

    storages_v = []
    scaled_renewables_v = Vector{Vector{Float64}}(undef,0)
    for (i, op) in enumerate(over_production)
        @info @sprintf("==== %d", i)
        @info @sprintf("storage_capacities = %s, over_production = %f", storage_capacities[i], op)
        scaled_renewables = get_scaled_renewables(powers, op)
        storages, scaled_renewables = compute_storage_level(powers.Load, scaled_renewables, storage_capacities[i], op, SF1_factor)
        push!(storages_v, storages)
        push!(scaled_renewables_v, scaled_renewables)
    end
    storages_v, scaled_renewables_v
end

"""
    determine_overproduction(Load, RP)

    determine minimum storage capacity as a function of over production

    Load : detrended Load
    RP : detrended and scaled renewable power
    punit : unit of Load and RP (MW. GW, TW)
"""
function determine_overproduction(dates::Vector{DateTime}, Load::Vector{Float64}, WWSB_sc::Vector{Float64}, punit::String)
    overproduction = collect(LinRange(1.05, 1.5, 20))
    storage_capacities = []
    for op in overproduction
        stc1 = 1.0
        stc2 = 1.0
        minS = -1.0
        it = 0
        while minS < 0.0 && it < 50

            storages, WWSB_sc = compute_storage_level(Load, WWSB_sc, storage_capacities, op::Float64, SF1_factor; log_p=false)

            min_storage_level1 = minimum(storages[1].SF)
            min_storage_level2 = minimum(storages[2].SF)

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
    compute_and_plot(par, power_data, power_data_de, storage_capacities, over_production)

"""
function compute_and_plot(par, power_data, power_data_de, storage_capacities, over_production)

    storages_v, WWSB_scaled = compute_storage_fill_level(power_data_de, storage_capacities, over_production, par.SF1_factor)
    
    dates = power_data.dates

    nb_stg = length(storages_v[1])
    @info "nb_stg =", nb_stg
    stores = []
    for j in 1:nb_stg
        push!(stores, [])
    end
    for j in 1:nb_stg
        for storages in storages_v
            push!(stores[j], storages[j])
        end
    end

    if par.plot_p
        @info @sprintf("fig_dir = %s", par.fig_dir)

        fig      = [1]
        Load     = power_data.Load
        Load_de  = power_data_de.Load
        WWSB_de  = power_data_de.WWSBPower

        ΔEL = (WWSB_de - Load_de)
        
        plot_powers(dates, Load, Load_de, power_data.WWSBPower, WWSB_de, 0, par.fig_dir, par.punit, fig)

        plot_detrended(dates, power_data.WWSBPower, WWSB_de, ΔEL, 
            Load, Load_de, power_data_de.Load_trend, par.punit, par.fig_dir, fig, data_are_averaged = false)
        
        plot_cumulative_power(dates, WWSB_de, Load_de, over_production, par.punit, par.fig_dir, fig)

        for j in 1:nb_stg
            @info j, fig
            plot_storage_fill_level(dates, Load_de, WWSB_de, WWSB_scaled, stores[j], 
                over_production, j, par.fig_dir, fig, par.punit, plot_all_p = par.plot_all_p)
        end
    end
end

"""
    load data and compute and plot storage fill levels, data are smoothed using moving averages
"""
function compute_and_plot_averaged(storage_capacities::Vector{Float64}, over_production::Vector{Float64}, par)

    power_data = PowerData(par);
    power_data_av = get_averaged_power_data(power_data, par.averaging_hours, par.averaging_method)

    power_data_avde = if par.scale_to_installed_power_p
        renewables_detrend_and_scale(par, power_data_av);
    else
        detrend_renewables(power_data_av);
    end

    storages_v, WWSB_scaled = compute_storage_fill_level(power_data_avde, storage_capacities, over_production, par.SF1_factor)
    nb_stg = length(storages_v[1])

    nb_stg = length(storages_v[1])
    @info @sprintf("nb_storages = %s", nb_stg)
    stores = []
    for j in 1:nb_stg
        push!(stores, [])
    end
    for j in 1:nb_stg
        for storages in storages_v
            push!(stores[j], storages[j])
        end
    end

    if par.plot_p
        @info @sprintf("fig_dir = %s", par.fig_dir)

        fig = [1]
        Load    = power_data_av.Load
        WWSB    = power_data_av.WWSBPower
        Load_de = power_data_avde.Load
        WWSB_de = power_data_avde.WWSBPower

        ΔEL = (WWSB_de - Load_de)
        dates = power_data_av.dates

        plot_powers(dates, Load, Load_de, WWSB, WWSB_de, 0, par.fig_dir, par.punit, fig)

        plot_detrended(dates, WWSB, WWSB_de, power_data_avde.WWSB_trend, ΔEL, Load, Load_de, 
                        power_data_avde.Load_trend, par.punit, par.fig_dir, fig, data_are_averaged = false)
                
        plot_cumulative_power(dates, WWSB_de, Load_de, over_production, par.punit, par.fig_dir, fig)
    
        for j in 1:nb_stg
            @info j, fig
            plot_storage_fill_level(dates, Load_de, WWSB_de, WWSB_scaled, stores[j], over_production, j, par.fig_dir, fig, par.punit, plot_all_p = par.plot_all_p)
        end
    end
end

#    for (stg1, stg2, op) in zip(res1, res2, oprod)
#        Load_per_year = powers_to_energy_per_year(dates, Load)
#        IF1  = powers_to_energy_per_year(dates, stg1.IF) / Load_per_year
#        OF1  = powers_to_energy_per_year(dates, stg1.OF) / Load_per_year
#        CT1  = powers_to_energy_per_year(dates, stg1.CT) / Load_per_year
#        IM1  = powers_to_energy_per_year(dates, stg1.IM) / Load_per_year
#        IF2  = powers_to_energy_per_year(dates, stg2.IF) / Load_per_year
#        OF2  = powers_to_energy_per_year(dates, stg2.OF) / Load_per_year
#        CT2  = powers_to_energy_per_year(dates, stg2.CT) / Load_per_year
#        IM2  = powers_to_energy_per_year(dates, stg2.IM) / Load_per_year
#        SF1  = stg1.SF[end] / Load_per_year
#        SF2  = stg2.SF[end] / Load_per_year
#
#        @info @sprintf("==== comp_and_plot ================================================================================")
#        @info @sprintf("op = %3.1f,  sc1 = %10.4e, sc2 = %10.4e, SF1[end] = %8.2e, SF2[end] = %8.2e, Load/y = %8.2e", op, stg1.SC, stg2.SC, stg1.SF[end], stg2.SF[end], Load_per_year)
#        @info @sprintf("IF1 = %10.4e, OF1 = %10.4e, IF1-OF1 = %10.4e, SF1[end] = %10.4e, CT1 = %10.4e, IM1 = %10.4e", IF1, OF1,  IF1-OF1, SF1, CT1, IM1)
#        @info @sprintf("IF2 = %10.4e, OF2 = %10.4e, IF2-OF2 = %10.4e, SF2[end] = %10.4e, CT2 = %10.4e, IM2 = %10.4e", IF2, OF2,  IF1-OF2, SF2, CT2, IM2)
#    end
#
#    if plot_p
#        fig_dir = joinpath(fig_dir, "averaged")
#        mkpath(fig_dir)
#        fig = [1]
#        plot_powers(dates, Load_ec, Load_av, RP_ec, RP_av, averaging_hours, fig_dir, punit, fig)
#        plot_detrended(dates, RP, RP_av_de, RP_av_trend, ΔEL_av, Load_av, Load_av_de, Load_av_trend, punit, fig_dir, fig, data_are_averaged = true)
#        plot_storage_fill_level(dates, Load_av_de, RP_av_de, res1, oprod, "S1", fig_dir, fig, "storage_fill1", punit, plot_all_p = plot_all_p)
#        plot_storage_fill_level(dates, Load_av_de, RP_av_de, res2, oprod, "S2", fig_dir, fig, "storage_fill2", punit, plot_all_p = plot_all_p)
#    end
#end
