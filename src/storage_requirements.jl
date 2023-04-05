using Common

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("storage_requirements_functions.jl")

root = dirname(@__DIR__)
get_data_dir()   = joinpath(root, "data")
get_fig_dir()    = joinpath(root, "figures")
mkpath(get_fig_dir())

#"python.terminal.activateEnvironment": false

function get_storage_capacities(punit)

    stc2  = [14.0, 26.0, 35.0, 45.0]
    oprod = [1.5, 1.2, 1.15, 1.1, 1.05]

    stc2  = [1.0]
    oprod = [3.0]

    # stc1 < stc2
    stc1 = @. stc2 * 0.01
    stc2 = @. stc2 - stc1

    conversion_factor = uconvert("TW", punit)# * 0.6
    stc1 = stc1 .* conversion_factor
    stc2 = stc2 .* conversion_factor

    stc1, stc2, oprod
end

punit = ["MW", "GW", "TW"][2]
stc1, stc2, oprod = get_storage_capacities(punit)

plot_p     = [false, true][2]
plot_all_p = [false, true][1]
do_log     = [false, true][1]

comp_and_plot(stc1, stc2, oprod, get_data_dir(), get_fig_dir(), punit, plot_p=plot_p, plot_all_p=plot_all_p, do_log=do_log);
