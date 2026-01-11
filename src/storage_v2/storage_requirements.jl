include("../include_sisremo.jl")

import PyPlot as plt
plt.pygui(true)
plt.pygui(:qt5)

include("plot_results.jl")
include("storage_requirements_functions.jl")

#"python.terminal.activateEnvironment": false

"""
    get_storage_capacities(par, storage_caps)

    par - PowerParameter
    storage_caps - storage capacity
"""
function get_storage_capacities(par, storage_caps)
    storage_capacities = Vector{Vector{Float64}}(undef, 0)
    for sc in storage_caps
        # stc1 < stc2
        stc1 = @. sc * 0.01
        stc2 = @. sc
        if par.second_storage_p
            push!(storage_capacities, [stc1, stc2])
        else
            push!(storage_capacities, [stc2])
        end
    end
    storage_capacities
end

function make_parameter(start_year, end_year)
    par = PowerParameter()
    par.punit                        = u_GW # [1u_MW, 1u_GW, 1u_TW][3]
    par.scale_Bio                    = 1.0
    par.SF1_factor                   = 0.5
    par.scale_to_installed_power_p   = true
    par.plot_p                       = [false, true][2]
    par.plot_all_p                   = [false, true][1]
    par.log_p                        = [false, true][1]
    par.scale_with_installed_power_p = true
    par.averaging_hours              = 24*7*4
    par.averaging_method             = [:moving_average, :mean][1]
    par.start_year                   = start_year
    par.end_year                     = end_year
    par
end

function storage_and_overproduction(par)
    factor = uconversion_factor(par.punit, 1u_TW)
    storage_capacities = [x*factor for x in [14.0, 26.0, 35.0, 45.0, 55.0]]
    over_production = [1.5, 1.2, 1.15, 1.1, 1.05]
    storage_capacities, over_production
end

par = make_parameter(2016, 2025)

storage_capacities, over_production = storage_and_overproduction(par)

date1 = DateTime("2017-01-01")
date2 = DateTime("2025-12-31")

par = PowerParameter()
par.scale_with_installed_power_p = true

public_power = get_public_public_power(date1, date2, par)
plt.plot(public_power.uts, public_power.Won .* public_power.Woff)
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


plt.plot(public_power.Load)
plt.plot(detrended_power.Load)
plt.plot(detrended_power.Load_trend)

compute_and_plot(par, public_power, detrended_power, storage_capacities, over_production)

par.fig_dir = par.fig_dir*"_av4"
compute_and_plot_averaged(storage_capacities, over_production, par);
