import 

import PyPlot as pl
pl.pygui(true)
pl.pygui(:qt5)

include("../power_parameter.jl")
include("../storage_requirements_functions_v2.jl")

#"python.terminal.activateEnvironment": false

"""
    get_storage_capacities(punit, stc)

    par - StorageParamete
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

function make_parameter()
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
    par
end

function storage_and_overproduction(punit)
    factor = uconversion_factor(punit, 1u_TW)
    storage_capacities = [x*factor for x in [14.0, 26.0, 35.0, 45.0, 55.0]]
    over_production = [1.5, 1.2, 1.15, 1.1, 1.05]
    storage_capacities, over_production
end

function get_paths()
    # project root directory
    sisremo_dir = dirname(dirname(@__DIR__))
    hdf5_dir = joinpath(sisremo_dir, "data")
    json_dir = joinpath(hdf5_dir, "json_downloads")

    mkpath(hdf5_dir)
    mkpath(json_dir)
    json_dir, hdf5_dir
end

function get_power_data()
    json_dir, hdf5_dir = get_paths()
    start_year = 2016
    end_year = 2025
    par = make_parameter()

    storage_capacities, over_production = storage_and_overproduction(par.punit)

    power_data = PowerData(hdf5_dir, start_year, end_year, par);
    detrended_and_scaled_data = renewables_detrend_and_scale(hdf5_dir, power_data, par);

    power_data_de = if par.scale_to_installed_power_p
        renewables_detrend_and_scale(hdf5_dir, power_data, par);
    else
        detrend_renewables(power_data);
    end

    #dates       :: Vector{DateTime} # 1
    #uts         :: Vector{Int64}    # 2
    #Load        :: Vector{Float64}  # 3
    #Woff        :: Vector{Float64}  # 4
    #Won         :: Vector{Float64}  # 5 
    #Solar       :: Vector{Float64}  # 6 
    #Bio         :: Vector{Float64}  # 7 
    #WWSBPower   :: Vector{Float64}  # 8 
    #Load_trend  :: Vector{Float64}  # 9
    #Woff_trend  :: Vector{Float64}  # 10 
    #Won_trend   :: Vector{Float64}  # 11 
    #Solar_trend :: Vector{Float64}  # 12 
    #Bio_trend   :: Vector{Float64}  # 13 
    #WWSB_trend  :: Vector{Float64}  # 14 

    power_data_de.WWSBPower
    power_data_de.dates

    power_data_de.uts[2] - power_data_de.uts[1]

    power_data_de, par
end
