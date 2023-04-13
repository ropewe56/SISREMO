using Common

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("storage_requirements_functions_v2.jl")

root = dirname(@__DIR__)
get_data_dir()   = joinpath(root, "data")
get_fig_dir()    = joinpath(root, "figures")
mkpath(get_fig_dir())

#"python.terminal.activateEnvironment": false

function get_storage_capacities(punit, stc, oprod)

    conversion_factor = uconvert("TW", punit)# * 0.6

    storage_capacities = Vector{Vector{Float64}}(undef, 0)
    for sc in stc
        # stc1 < stc2
        stc1 = @. sc * 0.01
        stc2 = @. sc - stc1
        stc1 = stc1 * conversion_factor
        stc2 = stc2 * conversion_factor
        #push!(storage_capacities, [stc1, stc2])
        push!(storage_capacities, [stc2])
    end

    storage_capacities, oprod
end

function make_parameter()
    punit = ["MW", "GW", "TW"][2]

    start_year, stop_year = 2016, 2022
    scale_Bio  = 1.0
    scale_to_installed_power_p = true

    plot_p     = [false, true][1]
    plot_all_p = [false, true][1]
    do_log     = [false, true][1]

    SF1_factor = 0.5

    par = Parameter(get_data_dir(), get_fig_dir(), punit, start_year, stop_year, scale_Bio, SF1_factor, scale_to_installed_power_p, plot_p, plot_all_p, do_log)

    stc   = [14.0, 26.0, 35.0, 45.0]
    oprod = [1.5, 1.2, 1.15, 1.1, 1.05]
    stc   = [1.0]
    oprod = [1.2]

    storage_capacities, oprod = get_storage_capacities(punit, stc, oprod)

    storage_capacities, oprod, par
end

storage_capacities, oprod, par = make_parameter()
comp_and_plot(storage_capacities, oprod, par);
