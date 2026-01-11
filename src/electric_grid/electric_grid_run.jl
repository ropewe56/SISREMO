import JLD
using Printf
using NLopt

import PyPlot as plt
plt.pygui(true)

include("../include_sisremo.jl")

#using Optimization
#using Optimization, OptimizationOptimJL
#using Optim
#import FiniteDifferences
#NoiseRobustDifferentiation.jl 

# /home/wester/Privat/Obsidian/Energiewende/Energie/Costs.md

#include("logger.jl")
include("electric_grid.jl")

#public_power, ppar = get_public_power();
#save_detrended_public_power(public_power, "save_detrended_public_power.hdf5")

const _TWh = 1.0e3

struct SimulationResult
    dates

    Load
    Prod
    Bat
    H2O
    Imp
    Exp
    Cur
    Res

    ELoad0 
    EProd0 

    sumELoad0
    sumEProd0
    sumLoadΔE
    sumProdΔEl
    sumProdΔEp
    sumBatΔEi 
    sumBatΔEo 
    sumH2OΔEi 
    sumH2OΔEo 
    sumImpΔE  
    sumExpΔE  
    sumCurΔE  
    sumResΔE  
    sumLoadC  
    sumProdCl 
    sumProdCp 
    sumBatCo  
    sumH2OCo  
    sumImpC   
    sumBatCi  
    sumH2OCi  
    sumExpC   
    sumCurC   
    sumResC   
end

function scale_values(years, Load0, Prod0, Load, Prod, Bat, H2O, Imp, Exp, Cur, Res)
    yearsTWh  = years * _TWh
    
    Load0    .= Load0    ./ yearsTWh
    Prod0    .= Prod0    ./ yearsTWh
    Load.ΔE  .= Load.ΔE  ./ yearsTWh
    Prod.ΔEl .= Prod.ΔEl ./ yearsTWh
    Prod.ΔEp .= Prod.ΔEp ./ yearsTWh
    Bat.ΔEi  .= Bat.ΔEi  ./ yearsTWh
    Bat.ΔEo  .= Bat.ΔEo  ./ yearsTWh
    H2O.ΔEi  .= H2O.ΔEi  ./ yearsTWh
    H2O.ΔEo  .= H2O.ΔEo  ./ yearsTWh
    Imp.ΔE   .= Imp.ΔE   ./ yearsTWh
    Exp.ΔE   .= Exp.ΔE   ./ yearsTWh
    Cur.ΔE   .= Cur.ΔE   ./ yearsTWh
    Res.ΔE   .= Res.ΔE   ./ yearsTWh

    Load.C   .= Load.C   ./ years
    Prod.Cl  .= Prod.Cl  ./ years
    Prod.Cp  .= Prod.Cp  ./ years
    Bat.Co   .= Bat.Co   ./ years
    H2O.Co   .= H2O.Co   ./ years
    Imp.C    .= Imp.C    ./ years
    Bat.Ci   .= Bat.Ci   ./ years
    H2O.Ci   .= H2O.Ci   ./ years
    Exp.C    .= Exp.C    ./ years
    Cur.C    .= Cur.C    ./ years
    Res.C    .= Res.C    ./ years
end

function SimulationResult(dates, ELoad0, EProd0, Load, Prod, Bat, H2O, Imp, Exp, Cur, Res)
    years = Dates.value(dates[end] - dates[1])/(3600*1.0e3)/(365.0*24.0)    
    #scale_values(years, ELoad0, EProd0, Load, Prod, Bat, H2O, Imp, Exp, Cur, Res)
    
    dy  = 1.0 / years
    dyW = 1.0 / (years*_TWh)

    sumELoad0 = sum(ELoad0) * dyW
    sumEProd0 = sum(EProd0) * dyW   

    sumLoadΔE  = sum(Load.ΔE) * dyW

    sumProdΔEl = sum(Prod.ΔEl) * dyW
    sumProdΔEp = sum(Prod.ΔEp) * dyW

    sumBatΔEi  = sum(Bat.ΔEi) * dyW
    sumBatΔEo  = sum(Bat.ΔEo) * dyW

    sumH2OΔEi  = sum(H2O.ΔEi) * dyW
    sumH2OΔEo  = sum(H2O.ΔEo) * dyW
    
    sumImpΔE   = sum(Imp.ΔE) * dyW
    sumExpΔE   = sum(Exp.ΔE) * dyW
    sumCurΔE   = sum(Cur.ΔE) * dyW
    sumResΔE   = sum(Res.ΔE) * dyW

    sumLoadC  = sum(Load.C)  * dy
    sumProdCl = sum(Prod.Cl) * dy
    sumProdCp = sum(Prod.Cp) * dy
    sumBatCo  = sum(Bat.Co)  * dy
    sumH2OCo  = sum(H2O.Co)  * dy
    sumImpC   = sum(Imp.C)   * dy

    sumBatCi  = sum(Bat.Ci) * dy
    sumH2OCi  = sum(H2O.Ci) * dy
    sumExpC   = sum(Exp.C)  * dy

    sumCurC   = sum(Cur.C) * dy
    sumResC   = sum(Res.C) * dy


    SimulationResult(dates, Load, Prod, Bat, H2O, Imp, Exp, Cur, Res,                    
        ELoad0, 
        EProd0, 
        sumELoad0 , 
        sumEProd0 , 
        sumLoadΔE , 
        sumProdΔEl, 
        sumProdΔEp, 
        sumBatΔEi , 
        sumBatΔEo , 
        sumH2OΔEi , 
        sumH2OΔEo , 
        sumImpΔE  , 
        sumExpΔE  , 
        sumCurΔE  , 
        sumResΔE  , 
        sumLoadC  , 
        sumProdCl , 
        sumProdCp , 
        sumBatCo  , 
        sumH2OCo  , 
        sumImpC   , 
        sumBatCi  , 
        sumH2OCi  , 
        sumExpC   , 
        sumCurC   , 
        sumResC   )
                
end

function print_results(sr)
    minH2O = minimum(sr.H2O.E)
    ΔH2O   = abs(sr.H2O.E[1] - sr.H2O.E[end])
    fval   = get_cent_kWh(sr)

    #dt = (t2-t1)*1.0e-9 + (t3-t2)*1.0e-9 + (t4-t3)*1.0e-9
    @info @sprintf("minH2O = %10.4e, ΔH2O = %10.4e, sumRes = %10.4e, sumCurt = %10.4e", minH2O, ΔH2O, sum(sr.Res.ΔE), sum(sr.Cur.ΔE))

    @printf("sumELoad0   = %11.4e\n", sr.sumELoad0)
    @printf("sumLoadΔE   = %11.4e\n", sr.sumLoadΔE)

    @printf("sumProd0    = %11.4e\n", sr.sumEProd0)
    @printf("sumProdΔEpl = %11.4e\n", sr.sumProdΔEl + sr.sumProdΔEp)
    @printf("sumProdΔEl  = %11.4e\n", sr.sumProdΔEl)
    @printf("sumProdΔEp  = %11.4e\n\n", sr.sumProdΔEp)

    ΔEBat = sr.Bat.E[end] - sr.Bat.E[1]
    ΔEH2O = sr.H2O.E[end] - sr.H2O.E[1]

    @printf("sumBatΔEo   = %11.4e\n", sr.sumBatΔEo )
    @printf("sumBatΔEi   = %11.4e\n", sr.sumBatΔEi )
    @printf("ΔEBat       = %11.4e\n\n", ΔEBat )

    @printf("sumH2OΔEo   = %11.4e\n", sr.sumH2OΔEo )
    @printf("sumH2OΔEi   = %11.4e\n", sr.sumH2OΔEi  )
    @printf("ΔEH2O       = %11.4e\n\n", ΔEH2O )
    
    @printf("sumImpΔE    = %11.4e\n",   sr.sumImpΔE )
    @printf("sumResΔE    = %11.4e\n\n", sr.sumResΔE )

    @printf("sumExpΔE    = %11.4e\n",   sr.sumExpΔE)
    @printf("sumCurΔE    = %11.4e\n",   sr.sumCurΔE)
    @printf("sumResΔE    = %11.4e\n\n", sr.sumResΔE)

    @printf("\n")

    @printf("sumLoadC    = %11.4e\n", sr.sumLoadC)

    @printf("sumProdCl   = %11.4e\n", sr.sumProdCl)
    @printf("sumProdCp   = %11.4e\n", sr.sumProdCp)
    @printf("sumBatCo    = %11.4e\n", sr.sumBatCo )
    @printf("sumH2OCo    = %11.4e\n", sr.sumH2OCo )
    @printf("sumImpC     = %11.4e\n", sr.sumImpC  )
    @printf("sumResC     = %11.4e\n", sr.sumResC  )

    @printf("sumBatCi    = %11.4e\n", sr.sumBatCi)
    @printf("sumH2OCi    = %11.4e\n", sr.sumH2OCi)
    @printf("sumExpC     = %11.4e\n", sr.sumExpC )
    @printf("sumCurC     = %11.4e\n", sr.sumCurC )

    @printf("\n")
    ϵ = 1.0e-12
    #@printf("c_load   = %11.4e\n", sr.C_Load / (sr.E_Load0+ϵ) * 1.0e-9)
    #@printf("c_prodl  = %11.4e\n", sr.C_Prodl/ (sr.E_Prodl+ϵ) * 1.0e-9)
    #@printf("c_prodp  = %11.4e\n", sr.C_Prodp/ (sr.E_Prodp+ϵ) * 1.0e-9)
    #@printf("c_bato   = %11.4e\n", sr.C_Bato / (sr.E_Bato +ϵ) * 1.0e-9)
    #@printf("c_H2o    = %11.4e\n", sr.C_H2Oo / (sr.E_H2Oo +ϵ) * 1.0e-9)
    #@printf("c_impp   = %11.4e\n", sr.C_Imp  / (sr.E_Imp  +ϵ) * 1.0e-9)
    #@printf("\n")
    #@printf("c_bati   = %11.4e\n", sr.C_Bati / (sr.E_Bati +ϵ) * 1.0e-9)
    #@printf("c_H2i    = %11.4e\n", sr.C_H2Oi / (sr.E_H2Oi +ϵ) * 1.0e-9)
    #@printf("c_expp   = %11.4e\n", sr.C_Exp  / (sr.E_Exp  +ϵ) * 1.0e-9)
    #@printf("c_curt   = %11.4e\n", sr.C_Cur  / (sr.E_Cur  +ϵ) * 1.0e-9)
    #@printf("c_res    = %11.4e\n", sr.C_Res  / (sr.E_Res  +ϵ) * 1.0e-9)
    #@printf("\n")

    # Load 
    Eload = sr.sumProdΔEl + sr.sumBatΔEo + sr.sumH2OΔEo + sr.sumImpΔE + sr.sumResΔE

    Esinks = sr.sumBatΔEi + sr.sumH2OΔEi + sr.sumExpΔE  + sr.sumCurΔE
    
    ΔEload = Eload - sr.sumLoadΔE

    ΔEsink = Esinks + sr.sumProdΔEl - sr.sumEProd0

    # Load cost
    Cload = sr.sumProdCl + sr.sumBatCo + sr.sumH2OCo + sr.sumImpC + sr.sumResC 
    # Production profit
    Csink = sr.sumBatCi + sr.sumH2OCi + sr.sumExpC + sr.sumCurC

    ΔCload = Cload - sr.sumLoadC
    ΔCsink = Csink - sr.sumProdCp

    @printf("ΔEload      = %11.4e\n", ΔEload)
    @printf("ΔEsink      = %11.4e\n", ΔEsink)
    @printf("ΔCload      = %11.4e\n", ΔCload)
    @printf("ΔCsink      = %11.4e\n", ΔCsink)
end

@inline function get_cent_kWh(sr)
    C1 = (sr.sumLoadC)/sr.sumELoad0 * 1.0e-9
    C2 = (sr.sumCurC)/sr.sumELoad0 * 1.0e-9
    C1 + C2
end

function plot_all(sr::SimulationResult)
    plt.figure()
    plt.plot(sr.dates,  sr.EProd0, label="Prod")
    plt.plot(sr.dates, -sr.ELoad0, label="Load")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates, sr.Bat.E, label="Bat")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates, sr.H2O.E, label="H2")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.Bat.ΔEi, label="E_Bat_i")
    plt.plot(sr.dates, -sr.Bat.ΔEo, label="E_Bat_o")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.H2O.ΔEi, label="E_H2_i")
    plt.plot(sr.dates, -sr.H2O.ΔEo, label="E_H2_o")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.Imp.ΔE, label="E_Imp")
    plt.plot(sr.dates, -sr.Exp.ΔE, label="E_Exp")
    plt.legend()
    sr = compute(x, public_power, nhours, p);
    abs(get_cent_kWh(sr))
    

    plt.figure()
    plt.plot(sr.dates, -sr.Cur.ΔE, label="E_Cur")
    plt.plot(sr.dates,  sr.Res.ΔE, label="E_Res")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.Prod.Cl, label="C_Prodl")
    plt.plot(sr.dates,  sr.Prod.Cp, label="C_Prodp")
    plt.plot(sr.dates, -sr.Load.C,  label="C_Load")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.Bat.Co, label="C_Bato")
    plt.plot(sr.dates, -sr.H2O.Co, label="C_H2Oo")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates, sr.Imp.C, label="C_Imp")
    plt.plot(sr.dates, sr.Exp.C, label="C_Exp")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.Cur.C, label="C_Cur")
    plt.plot(sr.dates, -sr.Res.C, label="C_Res")
    plt.legend()

end


Base.@kwdef mutable struct EnergyParameter{T}
    prod_cost_factor :: T = T(0.7)

    Pro_CMWh     :: T = T(60.0)
    
    Bat_CMWho    :: T = T(80.0)
    Bat_CMWhi    :: T = T(20.0)
    Bat_Pin      :: T = T(20.0)
    Bat_Pout     :: T = T(30.0)
    Bat_ηin      :: T = T(0.9)
    Bat_ηout     :: T = T(0.9)
    Bat_Einit    :: T = T(-1.0)

    Bat_C0o      :: T = T(10.0)
    Bat_C0i      :: T = T(1.0)

    H2O_CMWho    :: T = T(120.0)
    H2O_CMWhi    :: T = T(20.0)
    H2O_Pin      :: T = T(50.0)
    H2O_Pout     :: T = T(70.0)
    H2O_ηin      :: T = T(0.7)
    H2O_ηout     :: T = T(0.7)
    H2O_Einit    :: T = T(-1.0)

    H2O_C0o      :: T = T(80.0)
    H2O_C0i      :: T = T(10.0)

    Imp_CMWh     :: T = T(70.0)
    Imp_Pin      :: T = T(5.0)

    Exp_CMWh     :: T = T(5.0)
    Exp_Pout     :: T = T(10.0)

    Cur_CMWh     :: T = T(200.0)
    Res_CMWh     :: T = T(500.0)

    fcall        :: Int64 = 0
    gcall        :: Int64 = 0

    x            :: Vector{T} = zeros(T, 4)
    h2od         :: T = 0.0
    batd         :: T = 0.0
    minhso       :: T = 0.0
    MINEH2O      :: T = 100.0
    dnorm        :: T = 1.0e-4
end

function create_heatmap(public_power, nhours, p, bcap::T, hcap::T, op::T) where T
    cost = Matrix{T}(undef, length(bcap), length(hcap))
    for (i,bc) in enumerate(bcap)
        for (j,hc) in enumerate(hcap)
            p.Bat_Einit = bc
            p.H2_Einit = hc
            x = [bc, hc, op]
            cost[i,j] = compute(x, public_power, nhours, p)
            @printf("%d, %d, %f, %f, %e\n", i, j, bc, hc, cost[i,j])
        end
        #@printf("%d, %d, %14.8e\n", i, 100, cost[i,100])
    end
    JLD.save("cost.jld", "cost", cost)
    cost
end

function compute(x, public_power, nhours, p::EnergyParameter{T}) where T
    t1 = time_ns()

    bcap, hcap, op = x[1], x[2], x[3]

    if p.Bat_Einit < 0.0
        p.Bat_Einit = bcap*0.5
    end
    if p.H2O_Einit < 0.0
        p.H2O_Einit = hcap*0.5
    end
    
    load = Load(copy(public_power.Load));

    years = Dates.value(public_power.dates[end] - public_power.dates[1])/(3600*1.0e3)/(365.0*24.0)

    #price_of_MWh = (one(T) + (op-one(T)) * p.prod_cost_factor) * p.prod_CostGWh
    Eprod0 = public_power.WWSBPower .* op
    prod  = Production(copy(Eprod0), p.Pro_CMWh);

    #ΔP = @. prod0 - load.Et# + (p.H2O_Pout + p.Exp_Pout)
    #ΔPmin = minimum(ΔP)
    #p.H2O_Pout = -(ΔPmin + p.Exp_Pout)

    bat  = make_battery(nhours,  bcap, p.Bat_Pin, p.Bat_Pout, p.Bat_ηin, p.Bat_ηout, p.Bat_CMWhi, p.Bat_CMWho, p.Bat_C0i, p.Bat_C0o, p.Bat_Einit);
    H2   = make_hydrogen(nhours, hcap, p.H2O_Pin, p.H2O_Pout, p.H2O_ηin, p.H2O_ηout, p.H2O_CMWhi, p.H2O_CMWho, p.H2O_C0i, p.H2O_C0o, p.H2O_Einit);

    impp = Import(nhours, p.Imp_Pin, p.Imp_CMWh)
    expp = Export(nhours, p.Exp_Pout, p.Exp_CMWh)
    curt = Curtailment(nhours, p.Exp_CMWh)
    res  = ResidualLoad(nhours, p.Res_CMWh)
    
    Δt = one(T)

    t2 = time_ns()
    run_system(load, prod, bat, H2, impp, expp, curt, res, Δt)
    t3 = time_ns()
    sr = SimulationResult(public_power.dates, public_power.Load, Eprod0, load, prod, bat, H2, impp, expp, curt, res)
    t4 = time_ns()
    
    sr
end

T = Float64
arrowpath = joinpath(DATAROOT, "detrended_and_scaled_data.arrow")
public_power = load_from_arrow(arrowpath);
nhours = length(public_power[!,:dates])
p = EnergyParameter{Float64}();

x = [1.0e3, 1.5e4, 1.3]
sr = compute(x, DetrendedPowerData(public_power), nhours, p);
abs(get_cent_kWh(sr))

print_results(sr)

plot_all(sr)

plt.plot(cumsum(sr.Load.C) - cumsum(sr.Prod.Cl + sr.Bat.Co + sr.H2O.Co + sr.Imp.C + sr.Res.C))
plt.plot(sr.Load.C - (sr.Prod.Cp + sr.Bat.Co + sr.H2O.Co + sr.Imp.C + sr.Res.C))

a = @. ( sr.Load.C - (sr.Prod.Cl + sr.Bat.Co + sr.H2O.Co + sr.Imp.C + sr.Res.C))
i = argmax(a)
sr.Load.C[i]
sr.Prod.Cl[i]
sr.Prod.ΔEp[i]
sr.Bat.Co[i]
sr.H2O.Co[i]
sr.Imp.C[i]
sr.Res.C[i]
sum(sr.Cur.C)
sum(a)
plt.plot(a)
cumsum(sr.Load.C)

sum(sr.H2O.ΔEo)/10.0/1.0e3
plt.plot(sr.H2O.ΔEo)