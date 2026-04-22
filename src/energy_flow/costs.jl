struct GenerationCosts
    Solar_invest   :: Float64 # ct/GWh
    Wind_invest    :: Float64 # ct/GWh
    Storage_invest :: Float64 # ct/GWh
    Solar          :: Float64 # ct/GWh
    Wind           :: Float64 # ct/GWh
    Storage        :: Float64 # ct/GWh
end

function GenerationCosts()  
    ct_EUR = 1.0 # ct
    kW_GW  = 1.0e6 # kW
    year_hours = (365.0*24.0)

    # fix cost
    windinvest_GW = 1000.0 * ct_EUR * 1.0e6      # 1000.0 EUR/kWp * EUR/EUR * kWp/GW = EUR/GW
    lifetime      = 20.0 * year_hours            # 20.0 years * h/year = h
    Wind_invest   = windinvest_GW / lifetime     # EUR/GWh

    # 800.0 * EUR / kWp
    solinvest_GW = 800.0 * ct_EUR * 1.0e6        # 800.0 EUR/kWp * EUR/EUR * kWp/GW = EUR/GW
    lifetime     = 20.0 * year_hours             # 20.0 years * h/year = h
    Solar_invest = solinvest_GW / lifetime       # EUR/GWh

    storage_inv    = 200.0 * ct_EUR * 1.0e6        # 200 EUR/kWh * EUR/EUR * kWh/GWh  = EUR/GWh
    #Storage_invest = storage_inv / full_cycles    # EUR/GWh     

    Solar   = 5.0 / ct_EUR * kW_GW # ct/GWh
    Wind    = 6.0 / ct_EUR * kW_GW # ct/GWh
    Storage = 4.0 / ct_EUR * kW_GW # ct/GWh
    
    GenerationCosts(Solar_invest, Wind_invest, storage_inv, Solar, Wind, Storage)
end
const generation_cost = GenerationCosts()

op = 3.8615670025348665
storage_capacity = 2500.0

function determine_cost(op::Float64, energy_flow::EnergyFlow, powers::Powers, storage_capacity)
    nhours = Float64(length(energy_flow.toload))
    year_hours = (365.0*24.0)
    GW  = 1.0e6 # kW

    #L0 = sum(powers.Load) / nhours
    L = sum(energy_flow.toload)       / nhours # GW/h
    S = sum(energy_flow.tostore)      / nhours # GW/h
    F = sum(energy_flow.fromstorage)  / nhours # GW/h
    C = sum(energy_flow.tocurtail)    / nhours # GW/h

    power_tot = L + (F-S) + C # GW
    factor = 1.0 / (GW * power_tot)

    wind_installed  = (powers.max_Woff + powers.max_Won) * op      # GW
    wind_inv_cth    = wind_installed * generation_cost.Wind_invest # GW * EUR/GWh = EUR/h 
    #wind_inv_ctkWh  = wind_inv_cth * factor                       # EUR/kWh

    sol_installed   = powers.max_Solar * op                        # GW
    sol_inv_cth     = sol_installed * generation_cost.Wind_invest  # GW * EUR/GWh = EUR/h 
    #sol_inv_ctkWh  = sol_inv_cth * factor                         # EUR/kWh

    full_cycles     = 1500.0
    storage_cost    = storage_capacity * generation_cost.Storage_invest # GWh EUR/GWh = EUR
    storage_usage   = sum(energy_flow.tostore) / (full_cycles * storage_capacity) # GWh/h / GWh = 1/h
    storage_gen_cth = storage_cost * storage_usage

    wind_energy     = ((sum(powers.Woff) + sum(powers.Won)) * op + sum(powers.Bio)) # GJ
    wind_power      = wind_energy / nhours                         # GJ / h = GW
    wind_gen_cth    = wind_power * generation_cost.Wind            # GW * EUR/GWh = EUR/h
    #wind_gen_ctkWh  = wind_gen_cth * factor                       # EUR/kWh

    sol_energy      = sum(powers.Solar) * op   
    sol_energy_cost = sum(powers.Solar) * op * generation_cost.Solar # GWh * EUR/GWh = EUR  
    sol_power       = sol_energy / nhours                          # GJ / h = GW
    sol_gen_cth     = sol_power * generation_cost.Solar            # GW * EUR/GWh = EUR/h
    #sol_gen_ctkWh   = sol_gen_cth * factor

    cost_wwsb = (wind_inv_cth + wind_gen_cth + sol_inv_cth + sol_gen_cth) # EUR/h
    cost_wwsb / power_tot # EUR/GWh 100.0 / 1.0e6 

    storage_gen_cth
    
    total_cost[1] / power_tot ,

    (wind_inv_cth + wind_gen_cth + sol_inv_cth + sol_gen_cth)*factor, storage_gen_cth*factor
end
