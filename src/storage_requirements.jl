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

    #stc2  = [26.0]
    #oprod = [1.2]

    # stc1 < stc2
    stc1 = @. stc2 * 0.1
    stc2 = @. stc2 - stc1

    conversion_factor = uconvert("TW", punit) * 0.6
    stc1 = stc1 .* conversion_factor
    stc2 = stc2 .* conversion_factor

    stc1, stc2, oprod
end

punit = ["MW", "GW", "TW"][2]
stc1, stc2, oprod = get_storage_capacities(punit)
plot_p     = [false, true][2]
plot_all_p = [false, true][2]
do_log     = [false, true][1]
comp_and_plot(stc1, stc2, oprod, get_data_dir(), get_fig_dir(), punit, plot_p=plot_p, plot_all_p=plot_all_p, do_log=do_log);

# comp_and_plot_averaged(get_data_dir(), enunit, plot_p=true);
7.01260e+04/(364*26)

m = 252453600000.0
m/(365*24*60*60*1.0e3)

B = stg.IF[i] - stg.OF[i] - (stg.SF[i] - stg.SF[i-1])

S =  2.60e+03
ΔP = -1.00e+01, S =  2.59e+03, I =  0.00e+00, O = -1.00e+01, C =  0.00e+00, M =  0.00e+00, B =  2.01e+01