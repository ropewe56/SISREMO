const SISREMOROOT = dirname(dirname(@__DIR__))
const DATAROOT    = joinpath(SISREMOROOT, "data")
const FIGDIR      = joinpath(SISREMOROOT, "figures")
const JSONROOT    = joinpath(DATAROOT, "json")
const DBPATH      = joinpath(DATAROOT, "ise_data.sqlite")

using Dates
using Printf
using Optim
using DataFrames
using Arrow
using Statistics


import PyPlot as plt
plt.pygui(true)
plt.pygui(:qt5)

include("../init_logging.jl")
include("../energy_data/include_energy_data.jl")

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

function get_detrended_power(par::PowerParameter; load_arrow = false)
    if load_arrow == true
        path = joinpath(DATAROOT, "detrended_power.arrow")
        @info load_from_arrow, path
        return load_from_arrow(path)
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


mutable struct Storage
    ηin         :: Float64         # efficiency while charging
    ηout        :: Float64         # efficiency while discharging
    capacity    :: Float64         # storage capacity              [punit_h] TWh       SC
    fill_level  :: Vector{Float64} # storage fill level            [punit_h] TWh       SF
    toload      :: Vector{Float64} # direct flow to load           [punit] TW, GW, MW  I2
    tostorecurt :: Vector{Float64} # to storage and/or curtailment [punit] TW, GW, MW  I3
    tostore     :: Vector{Float64} # to storage                    [punit] TW, GW, MW  I4
    tocurtail   :: Vector{Float64} # to curtailment                [punit] TW, GW, MW  I5
    fromstorage :: Vector{Float64} # from storage to load          [punit] TW, GW, MW  I6
    other       :: Vector{Float64} # other sources to load         [punit] TW, GW, MW  I7
end

function Storage(SC::Float64, nb_steps::Int64; ηin = 1.0, ηout = 1.0)
    SF = zeros(Float64, nb_steps)
    I2 = zeros(Float64, nb_steps)
    I3 = zeros(Float64, nb_steps)
    I4 = zeros(Float64, nb_steps)
    I5 = zeros(Float64, nb_steps)
    I6 = zeros(Float64, nb_steps)
    I7 = zeros(Float64, nb_steps)
    Storage(ηin, ηout, SC, SF, I2, I3, I4, I5, I6, I7)
end

function power_step_(stg::Storage, L, P, istep)
    ΔSC = stg.SC - stg.SF[istep-1]
    ΔSC = stg.capacity - stg.fill_level[istep-1]
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
"""
    power_step(stg::Storage, L, P, i)

    L - Load
    P - power
"""
function power_step(storage::Storage, Load, WWSB, istep)
    L = Load[istep]
    P = WWSB[istep]

    ΔSC = storage.capacity - storage.fill_level[istep-1]
    k = 0
    if P - L > 0.0
        storage.toload[istep]      = L
        storage.tostorecurt[istep] = P - L
        storage.tostore[istep]     = min(storage.tostorecurt[istep], ΔSC/storage.ηin)
        storage.fill_level[istep]  = storage.fill_level[istep-1] + storage.tostore[istep]*storage.ηin
        storage.tocurtail[istep]   = max(0.0, storage.tostorecurt[istep]-storage.tostore[istep])  # curtailment
        k = 1
    else
        storage.toload[istep]      = P
        storage.fromstorage[istep] = min(storage.fill_level[istep-1] * storage.ηout, L - P)
        storage.fill_level[istep]  = storage.fill_level[istep-1] - storage.fromstorage[istep]/storage.ηout
        storage.other[istep]       = L - storage.toload[istep] - storage.fromstorage[istep]      # residual load
        k = 2
    end
    k
end


"""
    compute_storage_level(Load, WWSB, torage_capacity)
    Load = power.Load
"""
function compute_storage_level(Load, WWSB, storage_capacity)
    nb_steps  = length(WWSB)
    storage = Storage(storage_capacity, nb_steps)
    storage.fill_level[1] = storage_capacity

    istep = 2
    while istep <= nb_steps
        power_step(storage, Load, WWSB, istep)        
        istep += 1
    end
    #@info istep, storage.fill_level[istep]
    #plt.plot(storage.fill_level)
    
    storage
end

struct ScaledPower
    Load :: Vector{Float64}
    Woff :: Vector{Float64}
    Won  :: Vector{Float64}
    Solar:: Vector{Float64}
    Bio  :: Vector{Float64}

    mean_Load  :: Float64
    mean_Woff  :: Float64
    mean_Won   :: Float64
    mean_Solar :: Float64
    mean_Bio   :: Float64
    mean_WWS   :: Float64
end
function ScaledPower(Load, Woff, Won, Solar, Bio)
    mean_Load = mean(Load)
    mean_Woff = mean(Woff)
    mean_Won  = mean(Won)
    mean_Sol  = mean(Solar)
    mean_Bio  = mean(Bio)
    mean_WWS  = mean_Woff + mean_Won + mean_Sol
    ScaledPower(Load, Woff, Won, Solar, Bio, mean_Load, mean_Woff, mean_Won, mean_Sol, mean_Bio, mean_WWS)
end

function scale_wind_and_solar_by_load4(power::ScaledPower, x)
    WWS = @. power.Woff * x[2] + power.Won * x[3] + power.Solar * x[4]

    scale = (power.mean_Load - power.mean_Bio) / mean(WWS) * x[1]
    Woff  = power.Woff  .* scale
    Won   = power.Won   .* scale
    Solar = power.Solar .* scale

    @. Woff + Won + Solar + Bio
end

function scale_wind_and_solar_by_load1(power::ScaledPower, x)
    Woff  = power.Woff
    Won   = power.Won
    Solar = power.Solar
    Bio   = power.Bio

    # mean_L*op = mean_WWS * scale + mean_B
    scale = (power.mean_Load - power.mean_Bio) / power.mean_WWS * x[1]

    Woff  = power.Woff  .* scale
    Won   = power.Won   .* scale
    Solar = power.Solar .* scale

    @. (Woff + Won + Solar) + Bio
end


function get_WWSB_scaled(power::ScaledPower, storage_capacity, x)
    WWSB = if length(x) == 4
        scale_wind_and_solar_by_load4(power, x)
    else
        scale_wind_and_solar_by_load1(power, x)
    end
    WWSB
end

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function minimize_overproduction(power, storage_capacity; info = false)
    
    function ff(x)
        WWSB = get_WWSB_scaled(power, storage_capacity, x)
        Load = power.Load
        storage = compute_storage_level(Load, WWSB, storage_capacity)

        storage_min = minimum(storage.fill_level)
        @info x, storage_min, storage_min^2 - storage_capacity*0.1 + sum(storage.other)

        storage_min^2 - storage_capacity*0.1 + sum(storage.other)
    end

    x = [6.388]
    results = Optim.optimize(ff, x, NelderMead())
    
    op = Optim.minimizer(results)
    storage = compute_storage(power, storage_capacity, op)
    storage_min = minimum(storage.fill_level)
    if info
        @info storage_capacity, (storage_min * (1.0 - 0.1))^2 + sum(storage.other), Optim.minimizer(results), Optim.minimum(results)
    end

    storage_capacity, Optim.minimizer(results), Optim.minimum(results)
end

"""
    determine_overproduction(Load, RP)

    determine minimum storage capacity as a function of over production

    Load : detrended Load
    RP : detrended and scaled renewable power
    punit : unit of Load and RP (MW. GW, TW)
"""
function determine_overproduction(power::ScaledPower, storage_capacities; info = true)
    ops = []
    storage_capacity = storage_capacities[2]
    for storage_capacity in storage_capacities
        push!(ops, minimize_overproduction(power, storage_capacity; info = info))
    end
    ops
end

"""
    PowerParameter definded in energy_data/power_data.jl
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

function init_run(start_year=2017, end_year=2024)
    par = make_power_parameter(start_year, end_year)
    detrended_power = get_detrended_power(par; load_arrow = true)
    storage_capacities = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0] .* uconversion_factor(par.punit, 1u_TW)
    power = ScaledPower(detrended_power.Load, detrended_power.Woff, detrended_power.Won, detrended_power.Solar, detrended_power.Bio)
    power, storage_capacities
end

function run_simulation(start_year=2017, end_year=2024)

    determine_overproduction(power, storage_capacities)
    storage = compute_storage(power, storage_capacities[3], 3.0)
    storage_min = extrema(storage.SF)[1]

    plt.plot(storage.SF)
   
    (storage_min * (1.0 - 0.1))^2 + sum(storage.I7)
end

start_year, end_year = 2017, 2024
op = 2.0
st = compute_storage(Load, WWSB, op)

extrema(st.SF)