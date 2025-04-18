import JLD
using Printf
using NLopt
#using Optim
using BaseUtils

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
    Load0
    Prod0 
    Load
    Prod
    Bat
    H2O
    Imp
    Exp
    Cur
    Res

    dates

    E_Load0
    E_Load 
    E_Prod0
    E_Prod 
    E_Bato  
    E_H2Oo  
    E_Imp 
    E_Res

    E_Bati
    E_H2Oi 
    E_Exp 
    E_Cur 

    C_Load 
    C_Prodo
    C_Prodi
    C_Bato
    C_H2Oo
    C_Imp 
    C_Res

    C_Bati
    C_H2Oi
    C_Exp 
    C_Cur 
end

function SimulationResult(dates, Load0, Prod0, Load, Prod, Bat, H2O, Imp, Exp, Cur, Res)
    years = Dates.value(dates[end] - dates[1])/(3600*1.0e3)/(365.0*24.0)
    
    E_Load0 = sum(Load0)   / years/_TWh
    E_Load  = sum(Load.ΔE) / years/_TWh
    
    E_Prod0 = sum(Prod0)   /years/_TWh
    E_Prod  = sum(Prod.ΔE) /years/_TWh

    E_Bati  = sum(Bat.ΔEi) /years/_TWh
    E_Bato  = sum(Bat.ΔEo) /years/_TWh

    E_H2Oi  = sum(H2O.ΔEi)  /years/_TWh
    E_H2Oo  = sum(H2O.ΔEo)  /years/_TWh
    
    E_Imp   = sum(Imp.ΔE) /years/_TWh
    E_Exp   = sum(Exp.ΔE) /years/_TWh
    E_Cur   = sum(Cur.ΔE) /years/_TWh
    E_Res   = sum(Res.ΔE) /years/_TWh

    C_Load  = sum(Load.C) /years
    C_Prodo = sum(Prod.Co) /years
    C_Prodi = sum(Prod.Ci) /years
    C_Bato  = sum(Bat.Co) /years
    C_H2Oo  = sum(H2O.Co) /years
    C_Imp   = sum(Imp.C)  /years

    C_Bati  = sum(Bat.Ci) /years
    C_H2Oi  = sum(H2O.Ci) /years
    C_Exp   = sum(Exp.C)  /years

    C_Cur   = sum(Cur.C) /years
    C_Res   = sum(Res.C) /years


    SimulationResult(Load0, Prod0, 
                    Load, Prod, Bat, H2O, Imp, Exp, Cur, Res,
                    dates, 
                    E_Load0,
                    E_Load ,
                    E_Prod0,
                    E_Prod ,
                    E_Bato ,
                    E_H2Oo  ,
                    E_Imp ,
                    E_Res ,

                    E_Bati ,
                    E_H2Oi  ,
                    E_Exp ,
                    E_Cur ,

                    C_Load ,
                    C_Prodo ,
                    C_Prodi ,
                    C_Bato ,
                    C_H2Oo  ,
                    C_Imp ,
                    C_Res  ,

                    C_Bati ,
                    C_H2Oi  ,
                    C_Exp ,
                    C_Cur )
end

function print_results(sr)
    @printf("E_Load0  = %11.4e\n", sr.E_Load0)
    @printf("E_Load   = %11.4e\n", sr.E_Load )

    @printf("E_Prod0  = %11.4e\n", sr.E_Prod0)
    @printf("E_Prod   = %11.4e\n", sr.E_Prod )
    @printf("E_Bato   = %11.4e\n", sr.E_Bato )
    @printf("E_H2Oo   = %11.4e\n", sr.E_H2Oo  )
    @printf("E_Imp    = %11.4e\n", sr.E_Imp )
    @printf("E_Res    = %11.4e\n", sr.E_Res )

    @printf("E_Bati   = %11.4e\n", sr.E_Bati )
    @printf("E_H2Oi   = %11.4e\n", sr.E_H2Oi  )
    @printf("E_Exp    = %11.4e\n", sr.E_Exp )
    @printf("E_Cur    = %11.4e\n", sr.E_Cur )
    @printf("E_Res    = %11.4e\n", sr.E_Res )

    @printf("\n")

    @printf("C_Load   = %11.4e\n", sr.C_Load)

    @printf("C_Prodo  = %11.4e\n", sr.C_Prodo)
    @printf("C_Prodi  = %11.4e\n", sr.C_Prodi)
    @printf("C_Bato   = %11.4e\n", sr.C_Bato )
    @printf("C_H2Oo   = %11.4e\n", sr.C_H2Oo)
    @printf("C_Imp    = %11.4e\n", sr.C_Imp)
    @printf("C_Res    = %11.4e\n", sr.C_Res)

    @printf("C_Bati   = %11.4e\n", sr.C_Bati )
    @printf("C_H2Oi   = %11.4e\n", sr.C_H2Oi)
    @printf("C_Exp    = %11.4e\n", sr.C_Exp)
    @printf("C_Cur    = %11.4e\n", sr.C_Cur)

    @printf("\n")
    ϵ = 1.0e-12
    @printf("c_load   = %11.4e\n", sr.C_Load / (sr.E_Load0+ϵ) * 1.0e-9)
    @printf("c_prod   = %11.4e\n", sr.C_Prodo/ (sr.E_Prod +ϵ) * 1.0e-9)
    @printf("c_bato   = %11.4e\n", sr.C_Bato / (sr.E_Bato +ϵ) * 1.0e-9)
    @printf("c_H2o    = %11.4e\n", sr.C_H2Oo / (sr.E_H2Oo +ϵ) * 1.0e-9)
    @printf("c_impp   = %11.4e\n", sr.C_Imp  / (sr.E_Imp  +ϵ) * 1.0e-9)
    @printf("\n")
    @printf("c_bati   = %11.4e\n", sr.C_Bati / (sr.E_Bati +ϵ) * 1.0e-9)
    @printf("c_H2i    = %11.4e\n", sr.C_H2Oi / (sr.E_H2Oi +ϵ) * 1.0e-9)
    @printf("c_expp   = %11.4e\n", sr.C_Exp  / (sr.E_Exp  +ϵ) * 1.0e-9)
    @printf("c_curt   = %11.4e\n", sr.C_Cur  / (sr.E_Cur  +ϵ) * 1.0e-9)
    @printf("c_res    = %11.4e\n", sr.C_Res  / (sr.E_Res  +ϵ) * 1.0e-9)
    @printf("\n")

    # Load 
    Esources = sr.E_Prod + sr.E_Bato + sr.E_H2Oo + sr.E_Imp + sr.E_Res

    Esinks = sr.E_Bati + sr.E_H2Oi + sr.E_Exp  + sr.E_Cur 
    Esourcebalance = Esources - sr.E_Load

    Esinkbalance = Esinks + sr.E_Prod - sr.E_Prod0

    @printf("Energy source balance  = %11.4e\n", Esourcebalance)
    @printf("Energy sink   balance  = %11.4e\n", Esinkbalance)

    # Load cost
    Csources = sr.C_Prodo + sr.C_Bato + sr.C_H2Oo + sr.C_Imp + sr.C_Res 
    # Production profit
    Csinks   = sr.C_Prodi + sr.C_Bati + sr.C_H2Oi + sr.C_Exp + sr.C_Cur

    Cbalance = Csources - sr.C_Load
    @printf("Load cost balance  = %11.4e\n", Cbalance)
    @printf("Sink cost balance  = %11.4e\n", Csinks)

end

@inline function get_cent_kWh(sr)
    (sr.C_Load)/sr.E_Load0 * 1.0e-9
end

function plot_all(sr::SimulationResult)
    plt.figure()
    plt.plot(sr.dates,  sr.Prod0, label="Prod")
    plt.plot(sr.dates, -sr.Load0, label="Load")
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

    plt.figure()
    plt.plot(sr.dates, -sr.Cur.ΔE, label="E_Cur")
    plt.plot(sr.dates,  sr.Res.ΔE, label="E_Res")
    plt.legend()

    plt.figure()
    plt.plot(sr.dates,  sr.Prod.C, label="C_Prod")
    plt.plot(sr.dates, -sr.Load.C, label="C_Load")
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

    Cur_CMWh     :: T = T(500.0)
    Res_CMWh     :: T = T(300000.0)

    fcall        :: Int64 = 0
    gcall        :: Int64 = 0

    x            :: Vector{T} = zeros(T, 4)
    h2od         :: T = 0.0
    batd         :: T = 0.0
    minhso       :: T = 0.0
    MINEH2O      :: T = 100.0
    dnorm        :: T = 1.0e-4
end

function create_heatmap(power_data, nhours, p, bcap::T, hcap::T, op::T) where T
    cost = Matrix{T}(undef, length(bcap), length(hcap))
    for (i,bc) in enumerate(bcap)
        for (j,hc) in enumerate(hcap)
            p.Bat_Einit = bc
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

function compute(x, power_data, nhours, p::EnergyParameter{T}) where T
    t1 = time_ns()

    bcap, hcap, op = x[1], x[2], x[3]

    if p.Bat_Einit < 0.0
        p.Bat_Einit = bcap*0.5
    end
    if p.H2O_Einit < 0.0
        p.H2O_Einit = hcap*0.5
    end

    load = Load(copy(power_data.Load));

    years = Dates.value(power_data.dates[end] - power_data.dates[1])/(3600*1.0e3)/(365.0*24.0)

    #price_of_MWh = (one(T) + (op-one(T)) * p.prod_cost_factor) * p.prod_CostGWh
    prod0 = power_data.WWSBPower .* op
    prod  = Production(copy(prod0), p.Pro_CMWh);

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
    sr = SimulationResult(power_data.dates, power_data.Load, prod0, load, prod, bat, H2, impp, expp, curt, res)
    t4 = time_ns()
    
    minhso = minimum(sr.H2O.E) - 10.0
    h2d    = abs(sr.H2O.E[1] - sr.H2O.E[end]) - 10.0
    fval   = get_cent_kWh(sr)

    #dt = (t2-t1)*1.0e-9 + (t3-t2)*1.0e-9 + (t4-t3)*1.0e-9
    @infoe @sprintf("[%8.2f, %8.2f, %8.2f] %8.4f,  %10.4f,  %10.4f, %10.4f", x[1], x[2], x[3], fval, sum(sr.Res.ΔE), minhso, h2d)

    sr
end

function comp(x, power_data, nhours, p)
    sr = compute(x, power_data, nhours, p);

    minEH2O = minimum(sr.H2O.E)
    maxRes = maximum(sr.Res.ΔE)
    print_results(sr)
    
    @printf("\n")
    @printf("C kWh   = %6.4f\n", get_cent_kWh(sr))
    @printf("minEH2O = %8.2e\n", minEH2O)
    @printf("ΔH2O    = %8.2e\n", sr.H2O.E[1] - sr.H2O.E[end])
    @printf("maxRes  = %8.2e\n", maxRes)
    @printf("\n")

    sr
end

import FiniteDifferences
function make_funcs(power_data, nhours, p)print_results(sr)


    fdm = FiniteDifferences.central_fdm(3, 1)
    
    function fn(x, p)
        sr = compute(x, power_data, nhours, p)
        r = abs(get_cent_kWh(sr))
        @infoe x, r
        r
    end

    function fnd(x)
        sr = compute(x, power_data, nhours, p)
        r = abs(get_cent_kWh(sr))
        @infoe x, r
        r
    end

    function fngrad(G, x, p)
        δ = x .* 1.0e-6

        f0 = fn(x, p)
        x[1] = x[1] + δ[1]
        f1 = fn(x, p)
        x[1] = x[1] - δ[1]print_results(sr)

        x[2] = x[2] + δ[2]
        f2 = fn(x, p)
        x[2] = x[2] - δ[2]
        x[3] = x[3] + δ[3]
        f3 = fn(x, p)
        x[3] = x[3] - δ[3]

        G[1] = (f1-f0)/δ[1]
        G[2] = (f2-f0)/δ[2]
        G[3] = (f3-f0)/δ[3]

        @infoe x, G
        #G .= FiniteDifferences.grad(fdm, fnd, x)
        #@infoe x, G
    end

    function fncons(x, p)
        sr = compute(x, power_data, nhours, p)
        res[1] = abs(sr.H2O.E[1] - sr.H2O.E[end])
        #res[2] = minimum(sr.H2O.E)
    end
        
    fn, fngrad, fncons
end

function find_optimum(lb, ub, u0, power_data, nhours, p)

    opt = NLopt.Opt(:GN_AGS, 3)
    #opt = NLopt.Opt(:LN_BOBYQA, 4)

    NLopt.lower_bounds!(opt, lb)
    NLopt.upper_bounds!(opt, ub)

    objective_fn, constraint_fn1, constraint_fn2, constraint_fn3 = make_funcs(power_data, nhours, p)
    #NLopt.inequality_constraint!(opt, constraint_fn1, 1e-8)
    #NLopt.inequality_constraint!(opt, constraint_fn3, 1e-8)

    NLopt.xtol_rel!(opt, 1e-4)
    NLopt.min_objective!(opt, objective_fn)

    min_f, min_x, ret = NLopt.optimize(opt, u0)
    num_evals = NLopt.numevals(opt)

    sr = comp(min_x, power_data, nhours, p)

    @printf("objective value : %6.4f\n", min_f)
    @printf("solution        : %s   \n", min_x)
    @printf("solution status : %s   \n", ret)
    @printf("nb_func         : %d   \n", num_evals)

    min_x, min_f, sr
end

T = Float64
power_data = load_detrended_power_data("save_detrended_power_data.hdf5");
nhours = length(power_data.dates)
p = EnergyParameter{Float64}();


x = [2.0e3, 1.6e4, 1.4]
sr = compute(x, power_data, nhours, p);
abs(get_cent_kWh(sr))

print_results(sr)


min_x, min_f, sr = find_optimum(lb, ub, u0, power_data, nhours, p);

x = min_x .+ [0.0, 0.0, 0.2, 0.0]
sr = compute(x, power_data, nhours, p)

#plot_all(sr)


#using Optim
#inner_optimizer = NelderMead()
#sol1 = optimize(func, lb, ub, u0, Fminbox(inner_optimizer)) #f_reltol = 0.01))#, x_abstol=1.0e-5, iterations=200))
#u1 = Optim.minimizer(sol1)

using Optimization
using OptimizationNLopt

lb = [1.0e3]
ub = [4.0e4]
λ  = 0.5
u0 = [((1.0-λ)*lb[i] + λ*ub[i]) for i in 1:1]

function fn(u, p)
    x = [1.0e3, u[1], 1.5]
    sr = compute(x, power_data, nhours, p)
    r = abs(get_cent_kWh(sr))
    m = minimum(sr.H2O.E)
    @infoe x, r, m
    m
end

function fngrad(G, u, p)
    δ = u[1] * 1.0e-3
    
    f0 = fn(u, p)
    u[1] = u[1] + δ
    f1 = fn(u, p)
    u[1] = u[1] - δ

    G[1] = (f1-f0)/δ
    #G .= FiniteDifferences.grad(fdm, fnd, x)
    @infoe u, G
end

#fn, fngrad, fncons = make_funcs(power_data, nhours, p);
optf = OptimizationFunction(fn, grad=fngrad);#, cons=fncons);
prob = OptimizationProblem(optf, u0, p, lb = lb, ub = ub)
sol = solve(prob, Opt(:LD_LBFGS, 1))
G = zeros(3)
fngrad(G, u0, p)


plot_results(sr)

x = [1.0e3, 2.0e4, 1.25]
sr = compute(x, power_data, nhours, p)
r = abs(get_cent_kWh(sr))
G=zeros(1)
fngrad(G, [1.25], p)
u=zeros(1)
u[1] = 1.25


C_Load   =  3.4181e+10
C_Prod   =  2.9280e+10
C_Bato   =  3.4667e+09
C_H2Oo   =  1.1011e+09
C_Imp    =  1.2864e+09

C_Prod + 