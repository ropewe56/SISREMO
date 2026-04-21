struct GenerationCosts
    Solar_invest   :: Float64 # ct/GWh
    Wind_invest    :: Float64 # ct/GWh
    Storage_invest :: Float64 # ct/GWh
    Solar          :: Float64 # ct/GWh
    Wind           :: Float64 # ct/GWh
    Storage        :: Float64 # ct/GWh
end

function GenerationCosts()  
    EUR = 100.0 # ct
    GW  = 1.0e6 # kW
    year_hours = (365.0*24.0)

    solinvest_GW = 800.0 * EUR * GW              # ct/GW
    lifetime     = 20.0 * year_hours             # h
    Solar_invest = solinvest_GW / lifetime       # ct/GWh

    windinvest_GW = 1000.0 * EUR * GW            # ct/GW
    lifetime      = 20.0 * year_hours            # h
    Wind_invest   = windinvest_GW / lifetime     # ct/GWh

    cycles        = 1500.0
    stginv_GW     = 20.0 * EUR * GW             # ct/GW
    lifetime      = cycles / 365.0 * year_hours  # h   1 cicle per day
    Storage_invest= stginv_GW / lifetime         # ct/GWh

    Solar   = 5.0 * GW # ct/GWh
    Wind    = 6.0 * GW # ct/GWh
    Storage = 4.0 * GW # ct/GWh
    
    GenerationCosts(Solar_invest, Wind_invest, Storage_invest, Solar, Wind, Storage)
end
const generation_cost = GenerationCosts()

op = 3.8615670025348665
storage_capacity = 2500.0

function determine_cost(op::Float64, energy_flow::EnergyFlow, powers::Powers, storage_capacity)
    nhours = Float64(length(energy_flow.toload))

    #L0 = sum(powers.Load) / nhours
    L = sum(energy_flow.toload)       / nhours # GW/h
    S = sum(energy_flow.tostore)      / nhours # GW/h
    #F = sum(energy_flow.fromstorage)  / nhours # GW/h
    C = sum(energy_flow.tocurtail)    / nhours # GW/h

    power_tot = L + S + C # GW
    factor = 1.0 / (GW * power_tot)

    wind_installed  = (powers.max_Woff + powers.max_Won) * op      # GW
    wind_inv_cth    = wind_installed * generation_cost.Wind_invest # GW * ct/GWh = ct/h 
    wind_inv_ctkWh  = wind_inv_cth * factor                        # ct/kWh

    sol_installed   = powers.max_Solar * op                        # GW
    sol_inv_cth     = sol_installed * generation_cost.Wind_invest  # GW * ct/GWh = ct/h 
    sol_invest      = sol_inv_cth * factor                         # ct/kWh

    wind_energy     = ((sum(powers.Woff) + sum(powers.Won)) * op + sum(powers.Bio)) # GJ
    wind_power      = wind_energy / nhours                         # GJ / h = GW
    wind_gen_cth    = wind_power * generation_cost.Wind            # GW * ct/GWh = ct/h
    wind_gen_ctkWh  = wind_gen_cth * factor                        # ct/h

    sol_energy      = sum(powers.Solar) * op   
    sol_power       = sol_energy / nhours                          # GJ / h = GW
    sol_gen_cth     = sol_power * generation_cost.Solar            # GW * ct/GWh = ct/h
    sol_gen_ctkWh   = sol_gen_cth * factor

    storage_invest  = storage_capacity * generation_cost.Storage_invest * factor
    storage_gen     = sum(energy_flow.tostore) / nhours * generation_cost.Storage * factor

    (wind_inv_ctkWh + wind_gen_ctkWh + sol_invest + sol_gen_ctkWh) , (storage_invest + storage_gen)
end

