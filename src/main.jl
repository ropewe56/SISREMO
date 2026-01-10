include("include_sisremo.jl")
import PyPlot as plt
plt.pygui(true)

date1 = DateTime("2025-01-01")
date2 = DateTime("2025-12-31")
tables = ["total_power", "public_power", "cbpf", "installed_power"]
tables = ["public_power"]
#data = select_from_db(date1, date2, tables)

par = PowerParameter()

power_data = get_public_power_data(date1, date2, par)
save_to_arrow(power_data, "power_data_2025.arrow")
power_data = load_from_arrow("power_data_2025.arrow")
PowerData(power_data)

installed_power = get_installed_power_data(power_data, par);
save_to_arrow(installed_power, "installed_data_2025.arrow")
installed_power = load_from_arrow("installed_data_2025.arrow")
InstalledPowerData(installed_power)

averaging_hours, averaging_method = 24*30, :mean # :moving_average
averaged_power = get_averaged_power_data(power_data, averaging_hours, averaging_method)
