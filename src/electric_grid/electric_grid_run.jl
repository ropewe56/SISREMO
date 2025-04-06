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

function plot_all(power_data, WWSBPower, prod, bat, H2, impp, expp)
    plt.figure()
    plt.plot(power_data.dates, WWSBPower, label="WWSB")
    plt.plot(power_data.dates, power_data.Load, label="L")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, bat.E, label="Bat")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, H2.E, label="H2")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, bat.ΔEi, label="Bi")
    plt.plot(power_data.dates, -bat.ΔEo, label="Bo")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, H2.ΔEi, label="Hi")
    plt.plot(power_data.dates, -H2.ΔEo, label="Ho")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, impp.ΔE, label="Import")
    plt.plot(power_data.dates, -expp.ΔE, label="Export")
    plt.legend()

    plt.figure()
    plt.plot(power_data.dates, prod.C, label="pC")
    plt.legend()
    plt.figure()
    plt.plot(power_data.dates, bat.C, label="bC")
    plt.plot(power_data.dates, H2.C, label="h2C")
    plt.legend()
    plt.figure()
    plt.plot(power_data.dates, impp.C, label="impC")
    plt.legend()

end

Base.@kwdef mutable struct EnergyParameter{T}
    prod_cost_factor :: T = T(0.7)
    prod_CostGWh     :: T = T(60.0)
    bat_CostGWh      :: T = T(80.0)
    H2_CostGWh       :: T = T(200.0)
    import_CostGWh   :: T = T(70.0)
    export_CostGWh   :: T = T(-70.0)
    export_outP      :: T = T(5.0)
    bat_PowerIn      :: T = T(20.0)
    bat_PowerOut     :: T = T(30.0)
    H2_PowerIn       :: T = T(50.0)
    H2_PowerOut      :: T = T(50.0)
    bat_ηin          :: T = T(0.9)
    bat_ηout         :: T = T(0.9)
    H2_ηin           :: T = T(0.7)
    H2_ηout          :: T = T(0.7)
    bat_Einit        :: T = T(50.0)
    H2_Einit         :: T = T(1.0e4)
    fcall            :: Int64 = 0
    gcall            :: Int64 = 0
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

function compute(x, power_data, nhours, p::EnergyParameter{T}, ret=:all) where T
    bcap, hcap, op = x[1], x[2], x[3]

    p.bat_Einit = bcap*0.5
    p.H2_Einit = bcap*0.5

    load = Load(power_data.Load);

    #price_of_MWh = (one(T) + (op-one(T)) * p.prod_cost_factor) * p.prod_CostGWh
    WWSBPower = power_data.WWSBPower.*op
    prod = Production(WWSBPower, p.prod_CostGWh);

    bat  = make_battery(nhours,  bcap, p.bat_PowerIn, p.bat_PowerOut, p.bat_ηin, p.bat_ηout, p.bat_CostGWh, p.bat_Einit);
    H2   = make_hydrogen(nhours, hcap, p.H2_PowerIn, p.H2_PowerOut, p.H2_ηin, p.H2_ηout, p.H2_CostGWh, p.H2_Einit);

    impp = Import(nhours, p.import_CostGWh)
    expp = Export(nhours, p.export_outP, p.export_CostGWh)
    curt = Curtailment(nhours, p.export_CostGWh)
    
    Δt = one(T)
    run_system(load, prod, bat, H2, impp, expp, curt, Δt)

    total_cost = sum(prod.C)+sum(impp.C)+sum(bat.C)+sum(H2.C)+sum(expp.C)
    #println(@sprintf("(%8.2e  %8.2e  %8.2e)  %16.8e", x[1], x[2], x[3], load.total_cost))
    if ret == :all
        return abs(total_cost), WWSBPower, load, prod, bat, H2, impp, expp
    end
    abs(total_cost)
end

T = Float64
power_data = load_detrended_power_data("save_detrended_power_data.hdf5");
nhours = length(power_data.dates)
p = EnergyParameter{Float64}()

x = [1.5e3, 50.0e3, 1.1]
cost, WWSBPower, load, prod, bat, H2, impp, expp = compute(x, power_data, nhours, p);
println(cost)
sum(prod.C)
sum(bat.C)
sum(H2.C)
sum(impp.C)
sum(expp.C)

argmax(expp.ΔE)

cost2 = +sum(impp.C)+sum(bat.C)+sum(H2.C)
cost2+sum(expp.C)

plot_all(power_data, WWSBPower, prod, bat, H2, impp, expp)

prod.Et
#plt.plot(bat.ΔEi)

