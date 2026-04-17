include("../include_sisremo.jl")

import PyPlot as plt
plt.pygui(true)
plt.pygui(:qt5)

include("plot_results.jl")
include("storage_requirements_functions.jl")

#"python.terminal.activateEnvironment": false




par = make_parameter(2016, 2025)

storage_capacities, over_production = get_storage_and_overproduction(par)

date1 = DateTime("2017-01-01")
date2 = DateTime("2025-12-31")

par = PowerParameter()
par.scale_with_installed_power_p = true

public_power = get_public_power(date1, date2, par)
save_to_arrow(public_power, joinpath(DATAROOT, "public_power.arrow"))
#public_power = load_from_arrow("public_power.arrow")
#PowerData(public_power)

installed_power = get_installed_public_power(public_power, par);
save_to_arrow(installed_power, joinpath(DATAROOT, "installed_power.arrow"))
#installed_power = load_from_arrow("installed_power.arrow")
#InstalledPowerData(installed_power)

detrended_power = get_detrended_public_power(public_power, installed_power, par)
save_to_arrow(detrended_power, joinpath(DATAROOT, "detrended_power.arrow"))

averaging_hours, averaging_method = 24*30, :mean # :moving_average
averaged_power = get_averaged_public_power(detrended_power, averaging_hours, averaging_method)
save_to_arrow(averaged_power, joinpath(DATAROOT, "averaged_power.arrow"))


plt.figure()
plt.plot(public_power.uts, public_power.Won .* public_power.Woff)

plt.figure()
plt.plot(public_power.Load)
plt.plot(detrended_power.Load)
plt.plot(detrended_power.Load_trend)

compute_and_plot(par, public_power, detrended_power, storage_capacities, over_production)

par.fig_dir = par.fig_dir*"_av4"
compute_and_plot_averaged(par, public_power, detrended_power, averaged_power, storage_capacities, over_production)
