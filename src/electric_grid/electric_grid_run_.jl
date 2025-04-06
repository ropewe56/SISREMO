import PyPlot as plt
plt.pygui(true)

#using Optimization
using Optimization, OptimizationOptimJL
#using Optim

# /home/wester/Privat/Obsidian/Energiewende/Energie/Costs.md

include("electric_grid_ut.jl")

function compute(x, power_data, nhours, par::EnergyParameter{T}) where T
    bcap, hcap, op = x[1], x[2], x[3]

    load = Load(power_data.Load);

    price_of_MWh = (one(T) + (op-one(T)) * par.prod_cost_factor) * par.prod_CostGWh
    prod = Production(power_data.WWSBPower*op, op, price_of_MWh);
    
    bat  = make_battery(nhours,  bcap, par.bat_PowerIn, par.bat_PowerOut, par.bat_ηin, par.bat_ηout, par.bat_CostGWh, par.bat_Einit);
    H2   = make_hydrogen(nhours, hcap, par.H2_PowerIn, par.H2_PowerOut, par.H2_ηin, par.H2_ηout, par.H2_CostGWh, par.H2_Einit);

    impp = Import(nhours, par.import_CostGWh)
    expp = Export(nhours, par.export_CostGWh)
    
    Δt = one(T)
    run_system(load, prod, bat, H2, impp, expp, Δt)

    println(@sprintf("(%8.2e  %8.2e  %8.2e)  %16.8e", x[1], x[2], x[3], load.total_cost))
    
    abs(load.total_cost), load, prod, bat, H2, impp, expp
end

T = Float64
function get_fopt()
    power_data, ppar = get_power_data();
    nhours = length(power_data.dates)
    ep = EnergyParameter{T}()
    fu(x, ppar) = compute(x, power_data, nhours, ep)[1]
    fu, ep
end
fu, ep = get_fopt()

u0 = [T( 150.0), T(2.0e4), T(1.4)]
lb = [T(  10.0), T(2.0e3), T(1.0)]
ub = [T(1000.0), T(5.0e4), T(1.5)]

#sol = optimize(f, lb, ub, u0, Fminbox(NelderMead()), iterations=100) # Fminbox(NelderMead()))) #IPNewton())

problem = OptimizationProblem(OptimizationFunction(fu, Optimization.AutoForwardDiff()),  u0, ep;
                        lb = lb,
                        ub = ub,
                        lcons = nothing,
                        ucons = nothing,
                        sense = nothing)

sol = solve(problem, Optimization.LBFGS(), maxiters = 100000)#, maxtime = 1000.0)



x = Optim.minimizer(sol)

maxiters

x = [ 1.0, 2.0e4, 1.0]
f(x)
load.total_cost, load, prod, bat, H2, impp, expp = compute(x, power_data, nhours, par);
load.total_cost
prod.price_of_MWh

plt.plot(power_data.dates, power_data.WWSBPower, label="WWSB")
plt.plot(power_data.dates, power_data.Load, label="L")
plt.legend()

plt.figure()
plt.plot(power_data.dates, bat.E, label="Bat")
plt.legend()
plt.figure()
plt.plot(power_data.dates, H2.E, label="H2")
plt.legend()

plt.figure()
plt.plot(power_data.dates, impp.Et, label="Import")
plt.legend()
plt.figure()
plt.plot(power_data.dates, expp.Et, label="Export")
plt.legend()


plt.figure()
plt.plot(power_data.dates, power_data.Load, label="L")
plt.plot(power_data.dates, prod.ΔE, label="P")
plt.plot(power_data.dates, bat.ΔE, label="B")
plt.plot(power_data.dates, H2.ΔE, label="H")
plt.plot(power_data.dates, impp.ΔE, label="Import")
plt.plot(power_data.dates, expp.ΔE, label="Export")
plt.legend()

500*1.0e9 * 0.2