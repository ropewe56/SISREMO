using BaseUtils

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("storage_requirements_functions_v2.jl")

root = dirname(@__DIR__)
get_data_dir()   = joinpath(root, "data")
get_fig_dir()    = joinpath(root, "figures")
mkpath(get_fig_dir())

#"python.terminal.activateEnvironment": false

"""
    get_storage_capacities(punit, stc, over_production)

    punit - power unit
    stc - storage capacity in MWh
    over_production
"""
function get_storage_capacities(punit, stc, over_production)

    conversion_factor = uconversion_factor(punit, 1u_TW)

    storage_capacities = Vector{Vector{Float64}}(undef, 0)
    for sc in stc
        # stc1 < stc2
        stc1 = @. sc * 0.01
        stc2 = @. sc
        stc1 = stc1 * conversion_factor
        stc2 = stc2 * conversion_factor
        #push!(storage_capacities, [stc1, stc2])
        push!(storage_capacities, [stc2])
    end

    storage_capacities, over_production
end

function make_parameter(;start_year = 2016, stop_year = 2024)
    punit = u_GW # [1u_MW, 1u_GW, 1u_TW][3]

    scale_Bio  = 1.0
    scale_to_installed_power_p = true

    plot_p     = [false, true][2]
    plot_all_p = [false, true][1]
    do_log     = [false, true][1]

    SF1_factor = 0.5

    par = Parameter(get_data_dir(), get_fig_dir(), punit, start_year, stop_year, scale_Bio, SF1_factor, scale_to_installed_power_p, plot_p, plot_all_p, do_log)

    stc = [14.0, 26.0, 35.0, 45.0, 55.0]
    over_production = [1.5, 1.2, 1.15, 1.1, 1.05]

    #stc = [45.0]
    #over_production = [1.2]

    storage_capacities, over_production = get_storage_capacities(punit, stc, over_production)

    storage_capacities, over_production, par
end

storage_capacities, over_production, par = make_parameter(start_year = 2016, stop_year = 2024);
compute_and_plot(storage_capacities, over_production, par);
compute_and_plot_averaged(storage_capacities, over_production, par);
