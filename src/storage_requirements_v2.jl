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

function get_storage_capacities(punit)

    stc   = [14.0, 26.0, 35.0, 45.0]
    oprod = [1.5, 1.2, 1.15, 1.1, 1.05]

    stc   = [2.5]
    oprod = [2.0]

    conversion_factor = uconvert("TW", punit)# * 0.6

    storage_capacities = Vector{Vector{Float64}}(undef, 0)
    for sc in stc
        # stc1 < stc2
        stc1 = @. sc * 0.01
        stc2 = @. sc - stc1
        stc1 = stc1 * conversion_factor
        stc2 = stc2 * conversion_factor
        push!(storage_capacities, [stc1, stc2])
        #push!(storage_capacities, [stc2])
    end

    storage_capacities, oprod
end

punit = ["MW", "GW", "TW"][2]
storage_capacities, oprod = get_storage_capacities(punit)

start_year, stop_year = 2016, 2022

plot_p     = [false, true][2]
plot_all_p = [false, true][2]
do_log     = [false, true][1]

comp_and_plot(storage_capacities, oprod, get_data_dir(), get_fig_dir(), punit, start_year, stop_year, plot_p=plot_p, plot_all_p=plot_all_p, do_log=do_log);
