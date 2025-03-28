import PyPlot as plt
plt.pygui(true)

using Optimization
using Optim

include("electric_grid.jl")

T = Float64
function get_fopt()
    power_data, ppar = get_power_data();
    nhours = length(power_data.dates)
    ep = EnergyParameter{Float64}()
    f(x) = compute(x, power_data, nhours, ep)[1]
    f, ep
end
f, ep = get_fopt()

u0 = [ 150.0, 2.0e4, 1.4]
lb = [  10.0, 2.0e3, 1.0]
ub = [1000.0, 5.0e4, 1.5]

sol = optimize(f, lb, ub, u0, Fminbox(NelderMead()), iterations=100) # Fminbox(NelderMead()))) #IPNewton())

problem = OptimizationProblem(f,  u0, ep;
                        lb = lb,
                        ub = ub,
                        lcons = nothing,
                        ucons = nothing,
                        sense = nothing)

    x = Optim.minimizer(sol)

maxiters

sol = solve(problem, Optimization.LBFGS(), maxiters = 100000, maxtime = 1000.0)

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