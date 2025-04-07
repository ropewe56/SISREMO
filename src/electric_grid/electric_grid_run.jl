import PyPlot as plt
plt.pygui(true)

import JLD
using Printf

#using Optimization
#using Optimization, OptimizationOptimJL
#using Optim
#import FiniteDifferences
#NoiseRobustDifferentiation.jl 

# /home/wester/Privat/Obsidian/Energiewende/Energie/Costs.md

#include("logger.jl")
include("electric_grid.jl")

#power_data, ppar = get_power_data();
#save_detrended_power_data(power_data, "save_detrended_power_data.hdf5")

const _TWh = 1.0e3

function plot_all(power_data, WWSBPower, prod, bat, H2, impp, expp, curt)
    plt.figure()
    plt.plot(power_data.dates, WWSBPower, label="WWSB")
    plt.plot(power_data.dates, power_data.Load, label="Load")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, bat.E, label="Bat")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, H2.E, label="H2")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, bat.ΔEi, label="bat_i")
    plt.plot(power_data.dates, -bat.ΔEo, label="bat_o")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, H2.ΔEi, label="H2_i")
    plt.plot(power_data.dates, -H2.ΔEo, label="H2_o")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, impp.ΔE, label="Import")
    plt.plot(power_data.dates, -expp.ΔE, label="Export")
    plt.plot(power_data.dates, curt.ΔE, label="Curt")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, prod.C, label="prod_C")
    plt.legend()
    plt.figure()
    plt.plot(power_data.dates, bat.C, label="bat_C")
    plt.plot(power_data.dates, H2.C, label="H2_C")
    plt.legend()
    plt.figure()
    plt.plot(power_data.dates, impp.C, label="imp_C")
    plt.legend()

end

Base.@kwdef mutable struct EnergyParameter{T}
    prod_cost_factor :: T = T(0.7)

    prod_CostMWh     :: T = T(60.0)
    
    C_bati      :: T = T(80.0)
    C_H2i       :: T = T(200.0)
    C_bato      :: T = T(80.0)
    C_H2o       :: T = T(200.0)
    
    import_CostMWh   :: T = T(70.0)
    export_CostMWh   :: T = T(-50.0)
    curt_CostMWh     :: T = T(-50.0)

    export_outP      :: T = T(10.0)
    imp_inP          :: T = T(10.0)

    bat_PowerIn      :: T = T(20.0)
    bat_PowerOut     :: T = T(30.0)
    H2_PowerIn       :: T = T(50.0)
    H2_PowerOut      :: T = T(50.0)

    bat_ηin          :: T = T(0.9)
    bat_ηout         :: T = T(0.9)

    H2_ηin           :: T = T(0.7)
    H2_ηout          :: T = T(0.7)

    bat_Einit        :: T = T(-1.0)
    H2_Einit         :: T = T(-1.0)

    fcall            :: Int64 = 0
    gcall            :: Int64 = 0
end

struct SimulationResult
    C_load 
    C_tot  
    C_prod 
    C_bati
    C_bato
    C_H2i
    C_H2o
    C_impp 
    C_expp 
    C_curt 
    E_load0
    E_load 
    E_prod0
    E_prod 
    E_bati
    E_bato  
    E_H2i 
    E_H2o  
    E_impp 
    E_expp 
    E_curt 
end

function SimulationResult(power_data, WWSBPower, load, prod, bat, H2, impp, expp, curt)
    years = Dates.value(power_data.dates[end] - power_data.dates[1])/(3600*1.0e3)/(365.0*24.0)
    
    E_load0 = sum(power_data.Load)/_TWh/years
    E_load  = load.total_energy/_TWh/years

    C_load  = load.total_cost/years
    
    C_prod  = sum(prod.C)/years
    C_bato   = sum(bat.Co)/years
    C_H2o    = sum(H2.Co)/years
    C_bati   = sum(bat.Ci)/years
    C_H2i    = sum(H2.Ci)/years
    C_impp  = sum(impp.C)/years
    C_expp  = sum(expp.C)/years
    C_curt  = sum(curt.C)/years

    E_prod  = sum(prod.ΔE)/years/_TWh
    E_prod0 = sum(WWSBPower)/years/_TWh
    E_bati  = sum(bat.ΔEi) /years/_TWh
    E_bato  = sum(bat.ΔEo) /years/_TWh
    E_H2i   = sum(H2.ΔEi)  /years/_TWh
    E_H2o   = sum(H2.ΔEo)  /years/_TWh
    E_impp  = sum(impp.ΔE)/years/_TWh
    E_expp  = sum(expp.ΔE)/years/_TWh
    E_curt  = sum(curt.ΔE)/years/_TWh

    C_tot = C_prod + C_bato + C_H2o + C_impp#, C_expp, C_curt

    SimulationResult(C_load ,
                     C_tot  ,
                     C_prod ,
                     C_bati ,
                     C_bato ,
                     C_H2i  ,
                     C_H2o  ,
                     C_impp ,
                     C_expp ,
                     C_curt ,
                     E_load0,
                     E_load ,
                     E_prod0,
                     E_prod ,
                     E_bati ,
                     E_bato ,
                     E_H2i  ,
                     E_H2o  ,
                     E_impp ,
                     E_expp ,
                     E_curt )
end

@inline function get_cent_kWh(sr)
    sr.C_load/sr.E_load * 1.0e-9
end

function print_results(sr)
    @printf("E_load0  = %11.4e\n", sr.E_load0)
    @printf("E_load   = %11.4e\n", sr.E_load )
    @printf("E_prod0  = %11.4e\n", sr.E_prod0)
    @printf("E_prod   = %11.4e\n", sr.E_prod )
    @printf("E_bati   = %11.4e\n", sr.E_bati )
    @printf("E_bato   = %11.4e\n", sr.E_bato )
    @printf("E_H2i    = %11.4e\n", sr.E_H2i  )
    @printf("E_H2o    = %11.4e\n", sr.E_H2o  )
    @printf("E_impp   = %11.4e\n", sr.E_impp )
    @printf("E_expp   = %11.4e\n", sr.E_expp )
    @printf("E_curt   = %11.4e\n", sr.E_curt )
    @printf("\n")
    @printf("C_load   = %11.4e\n", sr.C_load)
    @printf("C_tot    = %11.4e\n", sr.C_tot )
    @printf("C_prod   = %11.4e\n", sr.C_prod)
    @printf("C_bati   = %11.4e\n", sr.C_bati )
    @printf("C_bato   = %11.4e\n", sr.C_bato )
    @printf("C_H2i    = %11.4e\n", sr.C_H2i)
    @printf("C_H2o    = %11.4e\n", sr.C_H2i)
    @printf("C_impp   = %11.4e\n", sr.C_impp)
    @printf("C_expp   = %11.4e\n", sr.C_expp)
    @printf("C_curt   = %11.4e\n", sr.C_curt)
    @printf("C_pbhiec = %11.4e\n", sr.C_prod + sr.C_bato + sr.C_H2o + sr.C_impp + sr.C_expp+sr.C_curt)
    @printf("\n")

    @printf("C_load   = %11.4e\n", sr.C_load/sr.E_load * 1.0e-9)
    @printf("C_prod   = %11.4e\n", sr.C_prod/sr.E_prod * 1.0e-9)
    @printf("C_bati   = %11.4e\n", sr.C_bati /sr.E_bati * 1.0e-9)
    @printf("C_bato   = %11.4e\n", sr.C_bato /sr.E_bato * 1.0e-9)
    @printf("C_H2i    = %11.4e\n", sr.C_H2i  /sr.E_H2i  * 1.0e-9)
    @printf("C_H2o    = %11.4e\n", sr.C_H2o  /sr.E_H2i  * 1.0e-9)
    @printf("C_impp   = %11.4e\n", sr.C_impp/sr.E_impp * 1.0e-9)
    @printf("C_expp   = %11.4e\n", sr.C_expp/sr.E_expp * 1.0e-9)
    #@printf("C_curt   = %11.4e\n", sr.C_curt/sr.E_curt)
end

function create_heatmap(power_data, nhours, p, bcap::T, hcap::T, op::T) where T
    cost = Matrix{T}(undef, length(bcap), length(hcap))
    for (i,bc) in enumerate(bcap)
        for (j,hc) in enumerate(hcap)
            p.bat_Einit = bc
            p.H2_Einit = hc
            x = [bc, hc, op]
            cost[i,j] = compute(x, power_data, nhours, p)
            @printf("%d, %d, %f, %f, %e\n", i, j, bc, hc, cost[i,j])
        end
        #@printf("%d, %d, %14.8e\n", i, 100, cost[i,100])
    end
    JLD.save("cost.jld", "cost", cost)
    cost
end

function compute(x, sr, nhours, p::EnergyParameter{T}, ret=:all) where T
    bcap, hcap, op = x[1], x[2], x[3]

    if p.bat_Einit < 0.0
        p.bat_Einit = bcap*0.5
    end
    if p.H2_Einit < 0.0
        p.H2_Einit = hcap*0.5
    end

    load = Load(power_data.Load);

    years = Dates.value(power_data.dates[end] - power_data.dates[1])/(3600*1.0e3)/(365.0*24.0)

    #price_of_MWh = (one(T) + (op-one(T)) * p.prod_cost_factor) * p.prod_CostGWh
    WWSBPower = power_data.WWSBPower.*op
    prod = Production(WWSBPower, p.prod_CostMWh);

    bat  = make_battery(nhours,  bcap, p.bat_PowerIn, p.bat_PowerOut, p.bat_ηin, p.bat_ηout, p.C_bati, p.C_bato, p.bat_Einit);
    H2   = make_hydrogen(nhours, hcap, p.H2_PowerIn,  p.H2_PowerOut,  p.H2_ηin,  p.H2_ηout,  p.C_H2i, p.C_H2o, p.H2_Einit);

    impp = Import(nhours, p.imp_inP, p.import_CostMWh)
    expp = Export(nhours, p.export_outP, p.export_CostMWh)
    curt = Curtailment(nhours, p.export_CostMWh)
    
    Δt = one(T)
    run_system(load, prod, bat, H2, impp, expp, curt, Δt)

    sr = SimulationResult(power_data, WWSBPower, load, prod, bat, H2, impp, expp, curt)
    
    if ret == :all
        return get_cent_kWh(sr), sr, load, prod, bat, H2, impp, expp, curt
    end

    abs(get_cent_kWh(sr))
end

T = Float64
power_data = load_detrended_power_data("save_detrended_power_data.hdf5");
nhours = length(power_data.dates)
p = EnergyParameter{Float64}()

x = [0.5_TWh, 30.0_TWh, 1.5]
cost, sr, load, prod, bat, H2, impp, expp, curt = compute(x, power_data, nhours, p);
print_results(sr)
get_cent_kWh(sr)

plot_all(power_data, sr.WWSBPower, prod, bat, H2, impp, expp, curt)

