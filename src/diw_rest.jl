#using JuMP
#using HiGHS
#model = Model(HiGHS.Optimizer)
#
#@variable(model, C2)
#@variable(model, C3)
#@variable(model, K1)
#
#@objective(model, Min, merit(x, v))
#
#@constraint(model, c1, constraint1(x1, v))
#@constraint(model, c2, constraint2(x2, v))
#@constraint(model, c3, constraint3(K1, v)

#ForwardDiff.gradient!(g, f, x)
#"""
#
#    Variables: Renewable capacity, Curtailment threshold, storage energy level[1]
#
#    Constraints
#        "$C$5"    = (sum(renewable_generation) + sum(discharging)) / storage_capacity) / 1000 * 100  = x1 Storage energy in GWh
#        "$C$6"    =  sum(curtailment)/sum(renewables_available)                        / 1000 * 100 <= x2 Renewable curtailment in per mille or percent
#        "$K$8777" = "$K$18"
#
#"""

x1 = [20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0] .* 1.0e-2;
x2 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
        110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500] ./ 1000.0;
x2