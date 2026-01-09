include("include_sisremo.jl")
import PyPlot as plt
plt.pygui(true)

date1 = DateTime("2025-01-01")
date2 = DateTime("2025-12-31")
tables = ["total_power", "public_power", "cbpf", "installed_power"]
tables = ["public_power"]
#data = select_from_db(date1, date2, tables)

par = PowerParameter()
power_data = PowerData(date1, date2, par);
installed_power = InstalledPowerData(power_data, par);

detrended_power = DetrendedPowerData(par, power_data, installed_power);

detrended_power_path = joinpath(DATAROOT, "detrended_power_2025.arrow")
save_detrended_power_data(detrended_power, detrended_power_path)

dp = load_detrended_power_data(detrended_power_path)
dp.dates

averaging_hours =  24*30
averaging_method = :moving_average
averaged_power = AveragedPowerData(detrended_power, averaging_hours, averaging_method)

plot_averaged(averaged_power, averaging_hours)