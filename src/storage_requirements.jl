using Common

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("storage_requirements_functions.jl")

root = dirname(@__DIR__)
get_data_dir()   = joinpath(root, "data")
get_fig_dir()    = joinpath(root, "figures")
mkpath(get_fig_dir())


eunit = ["MW", "GW", "TW"][2]

plot_p     = false#true
plot_all_p = false#true
comp_and_plot(get_data_dir(), get_fig_dir(), eunit, plot_p=plot_p, plot_all_p=plot_all_p);

# comp_and_plot_averaged(get_data_dir(), enunit, plot_p=true);
7.01260e+04/(364*26)

m = 252453600000.0
m/(365*24*60*60*1.0e3)
