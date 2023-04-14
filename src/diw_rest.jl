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
