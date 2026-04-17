

"""
    load data and compute and plot storage fille levels, original times (15 min)
    compute_and_plot(par, public_power, detrended_power, storage_capacities, over_production)

"""
function compute_and_plot(par, public_power, detrended_power, storage_capacities, over_production)

    storages_v, WWSB_scaled = compute_storage_fill_level(detrended_power, storage_capacities, over_production, par.SF1_factor)
    
    dates = public_power.dates

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
        fig_dir = par.fig_dir
        punit = par.punit 
        plot_all_p = par.plot_all_p

        fig     = [1]
        Load_ec = public_power.Load # energy charts,not detrended
        WWSB_ec = public_power.WWSBPower
        Load_de = detrended_power.Load
        WWSB_de = detrended_power.WWSBPower

        Load_trend = detrended_power.Load_trend

        ΔEL = (WWSB_de - Load_de)
        
        plot_powers(dates, Load_ec, Load_de, WWSB_ec, WWSB_de, 0, fig_dir, punit, fig)

        plot_detrended(dates, WWSB_ec, WWSB_de, ΔEL, 
            Load_ec, Load_de, Load_trend, punit, fig_dir, fig, data_are_averaged = false)
        
        plot_cumulative_power(dates, WWSB_de, Load_de, over_production, punit, fig_dir, fig)

        for j in 1:nb_stg
            @info j, fig
            plot_storage_fill_level(dates, Load_de, WWSB_de, WWSB_scaled, stores[j], 
                over_production, j, fig_dir, fig, punit, plot_all_p = plot_all_p)
        end
    end
end

"""
    load data and compute and plot storage fill levels, data are smoothed using moving averages
"""
function compute_and_plot_averaged(par::PowerParameter, 
        public_power::DataFrame, detrended_power::DataFrame, averaged_power::DataFrame, 
        storage_capacities, over_production);

    storages_v, WWSB_scaled = compute_storage_fill_level(averaged_power, storage_capacities, over_production, par.SF1_factor)
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
        punit = par.punit
        fig_dir = par.fig_dir

        fig = [1]
        Load_av    = Float64.(averaged_power.Load)
        WWSB_av    = Float64.(averaged_power.WWSBPower)
        dates_av   = DateTime.(averaged_power.dates)

        plot_averaged(dates_av, WWSB_av, Load_av, punit, fig_dir, fig)
                
        plot_cumulative_power(dates_av, WWSB_av, Load_av, over_production, par.punit, par.fig_dir, fig)
    
        for j in 1:nb_stg
            @info j, fig
            plot_storage_fill_level(dates, Load_av, WWSB_av, WWSB_scaled, stores[j], over_production, j, par.fig_dir, fig, par.punit, plot_all_p = par.plot_all_p)
        end
    end
end

