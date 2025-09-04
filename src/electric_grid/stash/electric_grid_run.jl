import JLD
using Printf

using SimpleLog

import PyPlot as plt
plt.pygui(true)

include("../energy_data/detrended_data.jl")

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

struct SimulationResult
    load0
    prod0 
    load
    prod
    bat
    H2
    impp
    expp
    curt
    res

    dates

    E_load0
    E_load 
    E_prod0
    E_prod 
    E_bato  
    E_H2o  
    E_impp 
    E_res

    E_bati
    E_H2i 
    E_expp 
    E_curt 

    C_load 
    C_prod 
    C_bato
    C_H2o
    C_impp 
    C_res

    C_bati
    C_H2i
    C_expp 
    C_curt 
end

function SimulationResult(dates, load0, prod0, load, prod, bat, H2, impp, expp, curt, res)
    years = Dates.value(dates[end] - dates[1])/(3600*1.0e3)/(365.0*24.0)
    
    E_load0 = sum(load0)   / years/_TWh
    E_load  = sum(load.ΔE) / years/_TWh
    
    E_prod0 = sum(prod0)   /years/_TWh
    E_prod  = sum(prod.ΔE) /years/_TWh

    E_bati  = sum(bat.ΔEi) /years/_TWh
    E_bato  = sum(bat.ΔEo) /years/_TWh

    E_H2i   = sum(H2.ΔEi)  /years/_TWh
    E_H2o   = sum(H2.ΔEo)  /years/_TWh
    
    E_impp  = sum(impp.ΔE) /years/_TWh
    E_expp  = sum(expp.ΔE) /years/_TWh
    E_curt  = sum(curt.ΔE) /years/_TWh
    E_res   = sum(res.ΔE)  /years/_TWh

    C_load  = sum(load.C) /years
    C_prod  = sum(prod.C) /years
    C_bato  = sum(bat.Co) /years
    C_H2o   = sum(H2.Co)  /years
    C_impp  = sum(impp.C) /years

    C_bati  = sum(bat.Ci) /years
    C_H2i   = sum(H2.Ci)  /years
    C_expp  = sum(expp.C) /years

    C_curt  = sum(curt.C) /years
    C_res   = sum(res.C) /years

    C_tot = C_prod + C_bato + C_H2o + C_impp#, C_expp, C_curt

    SimulationResult(load0, prod0, 
                    load, prod, bat, H2, impp, expp, curt, res,
                    dates, 
                    E_load0,
                    E_load ,
                    E_prod0,
                    E_prod ,
                    E_bato ,
                    E_H2o  ,
                    E_impp ,
                    E_res ,

                    E_bati ,
                    E_H2i  ,
                    E_expp ,
                    E_curt ,

                    C_load ,
                    C_prod ,
                    C_bato ,
                    C_H2o  ,
                    C_impp ,
                    C_res  ,

                    C_bati ,
                    C_H2i  ,
                    C_expp ,
                    C_curt )
end

function print_results(sr)
    @printf("E_load0  = %11.4e\n", sr.E_load0)
    @printf("E_load   = %11.4e\n", sr.E_load )

    @printf("E_prod0  = %11.4e\n", sr.E_prod0)
    @printf("E_prod   = %11.4e\n", sr.E_prod )
    @printf("E_bato   = %11.4e\n", sr.E_bato )
    @printf("E_H2o    = %11.4e\n", sr.E_H2o  )
    @printf("E_impp   = %11.4e\n", sr.E_impp )
    @printf("E_res    = %11.4e\n", sr.E_res )

    @printf("E_bati   = %11.4e\n", sr.E_bati )
    @printf("E_H2i    = %11.4e\n", sr.E_H2i  )
    @printf("E_expp   = %11.4e\n", sr.E_expp )
    @printf("E_curt   = %11.4e\n", sr.E_curt )
    @printf("E_res    = %11.4e\n", sr.E_res )

    @printf("\n")

    @printf("C_load   = %11.4e\n", sr.C_load)

    @printf("C_prod   = %11.4e\n", sr.C_prod)
    @printf("C_bato   = %11.4e\n", sr.C_bato )
    @printf("C_H2o    = %11.4e\n", sr.C_H2o)
    @printf("C_impp   = %11.4e\n", sr.C_impp)
    @printf("C_res    = %11.4e\n", sr.C_res)

    @printf("C_bati   = %11.4e\n", sr.C_bati )
    @printf("C_H2i    = %11.4e\n", sr.C_H2i)
    @printf("C_expp   = %11.4e\n", sr.C_expp)
    @printf("C_curt   = %11.4e\n", sr.C_curt)

    @printf("C_pbhiec = %11.4e\n", sr.C_prod + sr.C_bato + sr.C_H2o + sr.C_impp + sr.C_expp+sr.C_curt)
    @printf("\n")

    @printf("c_load   = %11.4e\n", sr.C_load / sr.E_load0* 1.0e-9)
    @printf("c_prod   = %11.4e\n", sr.C_prod / sr.E_prod * 1.0e-9)
    @printf("c_bato   = %11.4e\n", sr.C_bato / sr.E_bato * 1.0e-9)
    @printf("c_H2o    = %11.4e\n", sr.C_H2o  / sr.E_H2o  * 1.0e-9)
    @printf("c_impp   = %11.4e\n", sr.C_impp / sr.E_impp * 1.0e-9)
    @printf("c_bati   = %11.4e\n", sr.C_bati / sr.E_bati * 1.0e-9)
    @printf("c_H2i    = %11.4e\n", sr.C_H2i  / sr.E_H2i  * 1.0e-9)
    @printf("c_expp   = %11.4e\n", sr.C_expp / sr.E_expp * 1.0e-9)
    @printf("c_curt   = %11.4e\n", sr.C_curt / sr.E_curt * 1.0e-9)
    @printf("c_res    = %11.4e\n", sr.C_res  / sr.E_res  * 1.0e-9)
    @printf("\n")

    Esources = sr.E_bato + sr.E_H2o + sr.E_impp + sr.E_res   
    Esinks   = sr.E_bati + sr.E_H2i + sr.E_expp + sr.E_curt 
    Ebalance = sr.E_prod - Esinks + Esources - sr.E_load
    @printf("Ebalance  = %11.4e\n", Ebalance)

    Csources = sr.C_bato + sr.C_H2o + sr.C_impp + sr.C_res 
    Csinks   = sr.C_bati + sr.C_H2i + sr.C_expp + sr.C_curt
    Cbalance = Csources + Csinks + sr.C_prod - sr.C_load
    @printf("Cbalance  = %11.4e\n", Cbalance)

end

@inline function get_cent_kWh(sr)
    sr.C_load/sr.E_load0 * 1.0e-9
end

function plot_all(sr::SimulationResult)
    plt.figure()
    plt.plot(sr.dates, sr.prod0, label="Prod")
    plt.plot(sr.dates, sr.load0, label="Load")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates, sr.bat.E, label="Bat")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates, sr.H2.E, label="H2")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.bat.ΔEi, label="bat_i")
    plt.plot(sr.dates, -sr.bat.ΔEo, label="bat_o")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.H2.ΔEi, label="H2_i")
    plt.plot(sr.dates, -sr.H2.ΔEo, label="H2_o")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.impp.ΔE, label="Import")
    plt.plot(sr.dates, -sr.expp.ΔE, label="Export")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates, -sr.curt.ΔE, label="Curt")
    plt.plot(sr.dates,  sr.res.ΔE,  label="Residual")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.prod.C, label="prod_C")
    plt.plot(sr.dates, -sr.load.C, label="load_C")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates, sr.bat.Co, label="bat_C")
    plt.plot(sr.dates, sr.H2.Co, label="H2_C")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates, sr.impp.C, label="imp_C")
    plt.plot(sr.dates, sr.expp.C, label="exp_C")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates, -sr.curt.C, label="C_Curt")
    plt.plot(sr.dates,  sr.res.C,  label="C_Residual")
    plt.legend()

end


Base.@kwdef mutable struct EnergyParameter{T}
    prod_cost_factor :: T = T(0.7)

    prod_CostMWh     :: T = T(60.0)
    
    C_bato           :: T = T(80.0)
    C_bati           :: T = T(20.0)
    bat_PowerIn      :: T = T(20.0)
    bat_PowerOut     :: T = T(30.0)
    bat_Einit        :: T = T(-1.0)
    bat_ηin          :: T = T(0.9)
    bat_ηout         :: T = T(0.9)

    C_H2o            :: T = T(200.0)
    C_H2i            :: T = T(20.0)
    H2_PowerIn       :: T = T(50.0)
    H2_PowerOut      :: T = T(50.0)
    H2_ηin           :: T = T(0.7)
    H2_ηout          :: T = T(0.7)
    H2_Einit         :: T = T(-1.0)

    import_CostMWh   :: T = T(70.0)
    imp_inP          :: T = T(5.0)

    export_CostMWh   :: T = T(5.0)
    export_outP      :: T = T(10.0)

    curt_CostMWh     :: T = T(5.0)
    res_CostMWh      :: T = T(300.0)

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

    if p.bat_Einit < 0.0
        p.bat_Einit = bcap*0.5
    end
    if p.H2_Einit < 0.0
        p.H2_Einit = hcap*0.5
    end

    load = Load(copy(power_data.Load));

    years = Dates.value(power_data.dates[end] - power_data.dates[1])/(3600*1.0e3)/(365.0*24.0)

    #price_of_MWh = (one(T) + (op-one(T)) * p.prod_cost_factor) * p.prod_CostGWh
    prod0 = power_data.WWSBPower.*op
    prod  = Production(copy(prod0), p.prod_CostMWh);

    bat  = make_battery(nhours,  bcap, p.bat_PowerIn, p.bat_PowerOut, p.bat_ηin, p.bat_ηout, p.C_bati, p.C_bato, p.bat_Einit);
    H2   = make_hydrogen(nhours, hcap, p.H2_PowerIn,  p.H2_PowerOut,  p.H2_ηin,  p.H2_ηout,  p.C_H2i, p.C_H2o, p.H2_Einit);

    impp = Import(nhours, p.imp_inP, p.import_CostMWh)
    expp = Export(nhours, p.export_outP, p.export_CostMWh)
    curt = Curtailment(nhours, p.export_CostMWh)
    res  = ResidualLoad(nhours, p.res_CostMWh)
    
    Δt = one(T)
    run_system(load, prod, bat, H2, impp, expp, curt, res, Δt)

    sr = SimulationResult(power_data.dates, power_data.Load, prod0, load, prod, bat, H2, impp, expp, curt, res)
    
    if ret == :all
        return sr
    end

    ret = abs(get_cent_kWh(sr))
    println(ret, x)
    ret
end

T = Float64
power_data = load_detrended_power_data("save_detrended_power_data.hdf5");
nhours = length(power_data.dates)
p = EnergyParameter{Float64}()

x = [0.5_TWh, 30.0_TWh, 1.3]
x = u1
sr = compute(x, power_data, nhours, p);
print_results(sr)
get_cent_kWh(sr)

#plot_all(sr)

func(x) = compute(x, power_data, nhours, p, :single)


lb = [T(  10.0), T(2.0e3), T(1.3)]
ub = [T(1000.0), T(5.0e4), T(1.5)]

u0 = [0.5*(lb[i]+ub[i]) for i in 1:3]

using Optim
inner_optimizer = NelderMead()
sol1 = optimize(func, lb, ub, u0, Fminbox(inner_optimizer)) #f_reltol = 0.01))#, x_abstol=1.0e-5, iterations=200))
u1 = Optim.minimizer(sol1)
