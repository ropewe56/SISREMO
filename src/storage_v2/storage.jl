using Dates
using Printf

"""
    get the elapsed time of 1 step
    Dates.value(DateTime) => ms since 1 AD
    energy-chart data => 15 min time resolution ΔTh = 0.25
    1h = 60 * 60 * 1000 ms = 3.6e6 ms
    returns ΔTh - multiple of an hour per stpe
"""
function get_step_ΔTh(dates)
    ms = Dates.value(dates[2] - dates[1])
    ms/(3.6e6)
end

"""
    Get the number of years of the time series given by dates
"""
function number_of_years(dates)
    ΔT_ms_e  = Dates.value(dates[end] - dates[1])
    Δh_e     = ΔT_ms_e/(3.6e6)
    Δh_e/(365*24)
end

"""
    energy produced per year by power series P
"""
function powers_to_energy_per_year(dates, P)
    ΔTh = get_step_ΔTh(dates)
    nb_years = number_of_years(dates)
    sum(P)*ΔTh/nb_years
end


mutable struct Storage
    SC   :: Float64         # storage capacity              [punit_h] TWh
    SF   :: Vector{Float64} # storage fill level            [punit_h] TWh
    ηin  :: Float64         # efficiency while charging
    ηout :: Float64         # efficiency while discharging
    I2   :: Vector{Float64} # direct flow to load           [punit] TW, GW, MW
    I3   :: Vector{Float64} # to storage and/or curtailment [punit] TW, GW, MW
    I4   :: Vector{Float64} # to storage                    [punit] TW, GW, MW
    I5   :: Vector{Float64} # to curtailment                [punit] TW, GW, MW
    I6   :: Vector{Float64} # from storage to load          [punit] TW, GW, MW
    I7   :: Vector{Float64} # other sources to load         [punit] TW, GW, MW
end

function Storage(SC::Float64, nb_steps::Int64; ηin = 1.0, ηout = 1.0)
    SF = zeros(Float64, nb_steps)
    I2 = zeros(Float64, nb_steps)
    I3 = zeros(Float64, nb_steps)
    I4 = zeros(Float64, nb_steps)
    I5 = zeros(Float64, nb_steps)
    I6 = zeros(Float64, nb_steps)
    I7 = zeros(Float64, nb_steps)
    Storage(SC, SF, ηin, ηout, I2, I3, I4, I5, I6, I7)
end

function write_to_log(stg::Storage, ΔPin, out, i, j)
    B = stg.IF[i] - stg.I6[i] - (stg.SF[i] - stg.SF[i-1])
    write(out, @sprintf("%5d  %d  ΔP = %9.2e, S = %9.2e, I = %9.2e, O = %9.2e, C = %9.2e, M = %9.2e, B = %9.2e\n", i, j, ΔPin, stg.SF[i], stg.IF[i], stg.OF[i], stg.CT[i], stg.RE[i], B))
end
function write_power_step_to_log(stg::Storage, out, i, j)
    write(out, @sprintf("%5d  %d  P = %9.2e, L = %9.2e, I2 = %9.2e, I3 = %9.2e, I4 = %9.2e, I5 = %9.2e, I6 = %9.2e, I7 = %9.2e\n",
        i, j, stg.SF[i], P, L, stg.I2[i], stg.I3[i], stg.I4[i], stg.I5[i], stg.I6[i], stg.I7[i]))
end

"""
    power_step(stg::Storage, L, P, i)

    L - Load
    P - power
"""
function power_step(stg::Storage, L, P, istep)
    ΔSC = stg.SC - stg.SF[istep-1]
    k = 0
    if P - L > 0.0
        stg.I2[istep] = L
        stg.I3[istep] = P - L
        stg.I4[istep] = min(stg.I3[istep], ΔSC/stg.ηin)
        stg.SF[istep] = stg.SF[istep-1] + stg.I4[istep]*stg.ηin
        stg.I5[istep] = max(0.0, stg.I3[istep]-stg.I4[istep])  # curtailment
        k = 1
    else
        stg.I2[istep] = P
        stg.I6[istep] = min(stg.SF[istep-1] * stg.ηout, L - P)
        stg.SF[istep] = stg.SF[istep-1] - stg.I6[istep]/stg.ηout
        stg.I7[istep] = L - stg.I2[istep] - stg.I6[istep]      # residual load
        k = 2
    end
    k
end


function write_to_log(k, istep, j, load, scaled_renewables)
    if log_p
        out1 = open(joinpath(@__DIR__, "log1.log"), "w")
        out2 = open(joinpath(@__DIR__, "log2.log"), "w")

    write(out1, @sprintf("%3d, %d, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n", 
                istep, k, load, scaled_renewables,
                storages[j].I2[i],
                storages[j].I3[i],
                storages[j].I4[i],
                storages[j].I5[i],
                storages[j].I6[i],
                storages[j].I7[i],
                storages[j].SF[i],
                ))
    write(out2, @sprintf("%3d, %d, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n", 
                istep, k2, storages[j-1].I7[istep], storages[j-1].I5[istep],
                storages[j].I2[istep],
                storages[j].I3[istep],
                storages[j].I4[istep],
                storages[j].I5[istep],
                storages[j].I6[istep],
                storages[j].I7[istep],
                storages[j].SF[istep],
                ))
        
"""
    compute_storage_level(dates, Load, RP, punit, over_production, storage_capacity)

    given Load, RP, over_production and storage_capacity compute storage level as a funtion of time

    dates : times, ΔT = 1h
    Load  : power consumed
    RP    : renewable power production
    punit : power unit of Load and RP (MW, GW, TW)
    over_production  : renewable over production capacity factor, 1.0 is no over production capacity
    storage_capacity : storage capacity
"""
function compute_storage_level(load, scaled_renewables, storage_capacity, SF1_factor; log_p=false)
    nb_steps  = length(scaled_renewables)
    storage = Storage(storage_capacity, nb_steps)
    storage.SF[1] = storage_capacity*SF1_factor

    for istep in 2:nb_steps
        power_step(storage, load[istep], scaled_renewables[istep], istep)
    end
    
    storage
end

using Optim

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function compute_storage_fill_level(powers, storage_capacitiy, SF1_factor)
    function f(op)
        scaled_renewables = get_scaled_renewables(powers, op)
        storage = compute_storage_level(powers.Load, scaled_renewables, storage_capacitiy, SF1_factor)
        minimum(storage.SF)^2
    end
    op0 = [2.5]
    results = optimize(f, op0, NelderMead())
end

"""
    determine_overproduction(Load, RP)

    determine minimum storage capacity as a function of over production

    Load : detrended Load
    RP : detrended and scaled renewable power
    punit : unit of Load and RP (MW. GW, TW)
"""
function determine_overproduction(dates::Vector{DateTime}, Load::Vector{Float64}, WWSB::Vector{Float64}, punit::String)
    overproduction = collect(LinRange(1.05, 1.5, 20))
    storage_capacities = []
    for op in overproduction
        stc1 = 1.0
        stc2 = 1.0
        minS = -1.0
        it = 0
        while minS < 0.0 && it < 50

            storages, WWSB = compute_storage_level(Load, WWSB, storage_capacities, op::Float64, SF1_factor; log_p=false)

            min_storage_level1 = minimum(storages[1].SF)
            min_storage_level2 = minimum(storages[2].SF)

            stc1 = stc1 - min_storage_level1
            stc2 = stc2 - min_storage_level2
            it += 1
        end
        push!(storage_capacities, (stc1, stc2))
    end
end

"""
    get_overproduction_scaled_renewables(powers, op)

    powers : Load, Woff, Won, Solar, Bio
    op : over production factor

    R = Woff+Won+Solar
    s = (<Load>*op - <Bio>) / <Woff+Won+Solar>

    return (Woff+Won+Solar)/<Woff+Won+Solar> * (<Load>*op - <Bio>) + Bio
"""
function get_overproduction_scaled_renewables(powers, op)
    renewables = powers.Woff .+ powers.Won .+ powers.Solar
    mean_Load = mean(powers.Load)
    mean_Bio = mean(powers.Bio)
    mean_renewables = mean(renewables)

    # mean_L*op = mean_R * scale + mean_B
    scale = (mean_Load*op - mean_Bio) / mean_renewables

    renewables .* scale .+ powers.Bio
end

"""
    PoerParameter definded in energy_data/power_data.jl
"""
function make_power_parameter(start_year, end_year)
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

function get_detrended_power(par::PowerParameter; load_from_arrow = false)
    if load_from_arrow
        return load_from_arrow(joinpath(DATAROOT, "detrended_power.arrow"))
    end
    date1 = DateTime(par.start_year, 1, 1)
    date2 = DateTime(par.end_year, 12, 31)

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

    detrended_power
end

"""
    get_storage_capacities(par, storage_caps)

    par - PowerParameter
    storage_caps - storage capacity
"""
function get_storage_capacities(par, storage_caps)
    storage_capacities = Vector{Vector{Float64}}(undef, 0)
    for sc in storage_caps
        # stc1 < stc2
        stc1 = sc .* 0.01
        stc2 = sc
        if par.second_storage_p
            push!(storage_capacities, [stc1, stc2])
        else
            push!(storage_capacities, [stc2])
        end
    end
    storage_capacities
end

function get_storage_and_overproduction(par)
    factor = uconversion_factor(par.punit, 1u_TW)
    storage_capacities = [x*factor for x in [2.0, 3.0, 4.0, 8.0, 10.0]]
    over_production = [2.0, 4.0, 5.0]
    storage_capacities, over_production
end

function setup_simulation(start_year, end_year)

    par = make_power_parameter(start_year, end_year)
    detrended_power = get_detrended_power(par::PowerParameter; load_from_arrow = false)
    storage_capacities, over_production = get_storage_and_overproduction(par)
    
    storages_v, WWSB_scaled = compute_storage_fill_level(detrended_power, storage_capacities, over_production, par.SF1_factor)
