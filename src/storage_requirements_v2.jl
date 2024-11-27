import BaseUtils

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("storage_requirements_functions_v2.jl")

root = dirname(@__DIR__)
get_data_dir()   = joinpath(root, "data")
get_fig_dir()    = joinpath(root, "figures")
mkpath(get_data_dir())

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

function make_parameter(;start_year = 2016, end_year = 2024)

    par = StorageParameter()
    
    par.punit       = u_GW # [1u_MW, 1u_GW, 1u_TW][3]
    par.start_year  = start_year
    par.end_year   = end_year
    par.scale_Bio   = 1.0
    par.SF1_factor  = 0.5
    par.scale_to_installed_power_p = true
    par.second_storage_p = true
    par.plot_p      = [false, true][2]
    par.plot_all_p  = [false, true][1]
    par.log_p       = [false, true][1]
    par.scale_with_installed_power_p = true
    par.averaging_hours = 24*7*4

    factor = uconversion_factor(par.punit, 1u_TW)
    storage_capacities = [x*factor for x in [14.0, 26.0, 35.0, 45.0, 55.0]]
    over_production = [1.5, 1.2, 1.15, 1.1, 1.05]

    #stc = [45.0]
    #over_production = [1.2]

    storage_capacities = get_storage_capacities(par, storage_capacities)

    storage_capacities, over_production, par
end

storage_capacities, over_production, par = make_parameter(start_year = 2016, end_year = 2024);
compute_and_plot(storage_capacities, over_production, par);
par.fig_dir = par.fig_dir*"_av"
compute_and_plot_averaged(storage_capacities, over_production, par);
