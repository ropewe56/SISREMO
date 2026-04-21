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

mutable struct EnergyFlow
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

function EnergyFlow(capacity::Float64, nb_steps::Int64; ηin = 1.0, ηout = 1.0)
    fill_level  = zeros(Float64, nb_steps)
    toload      = zeros(Float64, nb_steps)
    tostorecurt = zeros(Float64, nb_steps)
    tostore     = zeros(Float64, nb_steps)
    tocurtail   = zeros(Float64, nb_steps)
    fromstorage = zeros(Float64, nb_steps)
    other       = zeros(Float64, nb_steps)
    EnergyFlow(ηin, ηout, capacity, fill_level, toload, tostorecurt, tostore, tocurtail, fromstorage, other)
end

"""
    power_step(stg::EnergyFlow, L, P, i)

    L - Load
    P - power
"""
function power_step(energy_flow::EnergyFlow, Load, WWSB, istep)
    L = Load[istep]
    P = WWSB[istep]

    ΔSC = energy_flow.capacity - energy_flow.fill_level[istep-1]
    k = 0
    if P - L > 0.0
        energy_flow.toload[istep]      = L
        energy_flow.tostorecurt[istep] = P - L
        energy_flow.tostore[istep]     = min(energy_flow.tostorecurt[istep], ΔSC/energy_flow.ηin)
        energy_flow.fill_level[istep]  = energy_flow.fill_level[istep-1] + energy_flow.tostore[istep]*energy_flow.ηin
        energy_flow.tocurtail[istep]   = max(0.0, energy_flow.tostorecurt[istep]-energy_flow.tostore[istep])  # curtailment
        k = 1
    else
        energy_flow.toload[istep]      = P
        energy_flow.fromstorage[istep] = min(energy_flow.fill_level[istep-1] * energy_flow.ηout, L - P)
        energy_flow.fill_level[istep]  = energy_flow.fill_level[istep-1] - energy_flow.fromstorage[istep]/energy_flow.ηout
        energy_flow.other[istep]       = L - energy_flow.toload[istep] - energy_flow.fromstorage[istep]      # residual load
        k = 2
    end
    k
end

"""
    compute_energy_flow(Load, WWSB, torage_capacity)
    Load = power.Load
"""
function compute_energy_flow(Load, WWSB_scaled, storage_capacity)
    nb_steps  = length(WWSB_scaled)
    energy_flow = EnergyFlow(storage_capacity, nb_steps)
    energy_flow.fill_level[1] = storage_capacity

    istep = 2
    while istep <= nb_steps
        power_step(energy_flow, Load, WWSB_scaled, istep)        
        istep += 1
    end
    
    energy_flow
end

