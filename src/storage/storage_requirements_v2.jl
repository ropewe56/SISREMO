import BaseUtils

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("power_parameter.jl")
include("storage_requirements_functions_v2.jl")

#"python.terminal.activateEnvironment": false

"""
    get_storage_capacities(punit, stc)

    par - StorageParamete
    storage_caps - storage capacity
"""
function get_storage_capacities(par, storage_caps)
    storage_capacities = Vector{Vector{Float64}}(undef, 0)
    for sc in storage_caps
        # stc1 < stc2
        stc1 = @. sc * 0.01
        stc2 = @. sc
        if par.second_storage_p
            push!(storage_capacities, [stc1, stc2])
        else
            push!(storage_capacities, [stc2])
        end
    end
    storage_capacities
end

function make_parameter()
    par = PowerParameter()
    par.punit                        = u_GW # [1u_MW, 1u_GW, 1u_TW][3]
    par.scale_Bio                    = 1.0
    par.SF1_factor                   = 0.5
    par.scale_to_installed_power_p   = true
    par.plot_p                       = [false, true][2]
    par.plot_all_p                   = [false, true][1]
    par.log_p                        = [false, true][1]
    par.scale_with_installed_power_p = true
    par.averaging_hours              = 24*7*4
    par.averaging_method             = [:moving_average, :mean][1]
    par
end

function storage_and_overproduction(par)
    factor = uconversion_factor(par.punit, 1u_TW)
    storage_capacities = [x*factor for x in [14.0, 26.0, 35.0, 45.0, 55.0]]
    over_production = [1.5, 1.2, 1.15, 1.1, 1.05]
    storage_capacities, over_production
end

function get_paths()
    # project root directory
    sisremo_dir = dirname(@__DIR__)
    hdf5_dir = joinpath(sisremo_dir, "data")
    json_dir = joinpath(hdf5_dir, "json_downloads")

    mkpath(hdf5_dir)
    mkpath(json_dir)
    json_dir, hdf5_dir
end
json_dir, hdf5_dir = get_paths()
start_year = 2016
end_year = 2025
par = make_parameter()

storage_capacities, over_production = storage_and_overproduction(par)

power_data = PowerData(hdf5_dir, start_year, end_year, par);
detrended_and_scaled_data = renewables_detrend_and_scale(hdf5_dir, power_data, par);

power_data_de = if par.scale_to_installed_power_p
    renewables_detrend_and_scale(hdf5_dir, power_data, par);
else
    detrend_renewables(power_data);
end


compute_and_plot(hdf5_dir, storage_capacities, over_production, par);

par.fig_dir = par.fig_dir*"_av4"
compute_and_plot_averaged(storage_capacities, over_production, par);
