include("include_sisremo.jl")
import PyPlot as plt
plt.pygui(true)

date1 = DateTime("2016-01-01")
date2 = DateTime("2025-12-31")

tables = ["total_power", "public_power", "cbpf", "installed_power"]
tables = ["public_power"]
#data = select_from_db(date1, date2, tables)

par = PowerParameter()
par.scale_with_installed_power_p = true

public_power = get_public_public_power(date1, date2, par)
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

