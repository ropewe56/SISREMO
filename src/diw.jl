using CSV
using DataFrames

root = "/home/wester/Projects/Julia/Private/Energy.jl/jupyter/Energy/DIW/eer-d-18-00186_code_data/model_and_data_spreadsheet_tool/"
csv = "spreadsheet_tool_1.csv"
data = CSV.read(joinpath(root, csv), DataFrame, header=4, delim=';', decimal=',');

mutable struct Vectors
    B :: Vector{Float64}
    C :: Vector{Float64}
    D :: Vector{Float64}
    E :: Vector{Float64}
    F :: Vector{Float64}
    G :: Vector{Float64}
    H :: Vector{Float64}
    I :: Vector{Float64}
    J :: Vector{Float64}
    K :: Vector{Float64}
end
function Vectors()
    B = Vector{Float64}(undef, 0)
    C = Vector{Float64}(undef, 0)
    D = Vector{Float64}(undef, 0)
    E = Vector{Float64}(undef, 0)
    F = Vector{Float64}(undef, 0)
    G = Vector{Float64}(undef, 0)
    H = Vector{Float64}(undef, 0)
    I = Vector{Float64}(undef, 0)
    J = Vector{Float64}(undef, 0)
    K = Vector{Float64}(undef, 0)
    Vectors(B, C, D, E, F, G, H, I, J, K)
end

ii = Dict(2012 => 3, 2013 => 5, 2014 => 7, 2015 => 9, 2016 => 11)
i = ii[2016]
v = Vectors()
v.B = data[:,i];
v.C = data[:,i+1];

"""

    Variables: Renewable capacity, Curtailment threshold, storage energy level[1]

    Constraints
        "$C$5"    = (sum(renewable_generation) + sum(discharging)) / storage_capacity) / 1000 * 100  = x1 Storage energy in GWh
        "$C$6"    =  sum(curtailment)/sum(renewables_available)                        / 1000 * 100 <= x2 Renewable curtailment in per mille or percent
        "$K$8777" = "$K$18"

"""
function minimize_max_storage_energy_level(v, x1, x2)
    #Variables: C2, C3, K1

    #minimize(maximum(K))
    #    C5 = (sum(F) + sum(J)) / sum(B)  = x1
    #    C6 =  sum(H) / sum(D)            <= x2
    #    K[end] = K1

    function constraint(v, x1, x2, K1)
        c1 = x1 - (sum(v.F) + sum(v.J)) / sum(v.B) # = 0
        c2 = x2 - sum(v.H) / sum(v.D)              # <= 0
        c1, c2, K1 - v.K[end]
    end


    OB = merit(v, C2, C3, K1)
    CS = constraint(v, x1, x2, K1)

end

# storage energy [GWh]
x1 = [20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0]

# Renewable curtailment [%]
x2 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
        110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500] ./ 100.0

maxK = [1.0e10, 0.0, 0.0, 0.0]

for SE in x1
    for CT in x2
        mK, sumL = minimize_max_storage_energy_level(B, C, SE, CT)
        if mK < maxK[1]
            maxK[1] = mK
            maxK[2] = sumL
            maxK[3] = SE
            maxK[4] = CT
        end
    end
end
println(maxK)


using JuMP
using HiGHS

function merit(v, C2, C3, K1)
    EfficiencyLoading     = 0.81
    EfficiencyDischarging = 0.926

    #B demand
    #C renewable availbility factor

    v.D = v.C .* C2                                 # renewable availability
    v.F = min.(v.D, v.B)                            # renewable generation
    v.E = v.D - v.F                                 # renewable surplus
    v.G = v.B - v.F                                 # residual demand
    v.H = min.(-v.E .+ C3, 0.0)                     # curtailment
    v.I = v.E - v.H;                                # storage loading

    n = length(v.D)
    v.J = zeros(Float64, n)                         # storage discharging
    v.J[1] = 0.0

    v.K = zeros(Float64, n)                         # storage energy level
    v.K[1] = K1

    for i in 2:n
        v.J[i] = min(v.G[i], v.K[i-1] .* EfficiencyDischarging);
        v.K[i] = v.K[i-1] + v.I[i] .* EfficiencyLoading - v.J[i] ./ EfficiencyDischarging;
    end

    v.L = map(x -> max(x, 0.0), v.G .- v.J)          # conventional generation
    v.M = v.B .- v.D;                                # residual demand (for RLDCs)

    maximum(v.K)
end

model = Model(HiGHS.Optimizer)

@variable(model, C2)
@variable(model, C3)
@variable(model, K1)

typeof(v)
@objective(model, Min, merit(v, C2, C3, K1))

@constraint(model, c1, x1 - (sum(v.F) + sum(v.J)) / sum(v.B) == 0.0)
@constraint(model, c2, x2 - sum(v.H) / sum(v.D) <=  0.0)
@constraint(model, c3, K1 == v.K[end])


using NLopt

function myfunc(x::Vector, grad::Vector)
    merit(v, x[1], x[2], x[3])
end

function myconstraint(x::Vector, grad::Vector, x1, x2)
    x1 - (sum(v.F) + sum(v.J)) / sum(v.B)
end

opt = Opt(:LD_LBFGS, 3)
opt.lower_bounds = [0.0, 0.0, 0.0]
opt.xtol_rel = 1e-4

const1(x1,v) = x1 - (sum(v.F) + sum(v.J)) / sum(v.B)
const2(x2,v) = x2 - sum(v.H) / sum(v.D)
const3(K1,v) =  K1 - v.K[end]

opt.min_objective = myfunc

equality_constraint!(opt, (x, g)   -> const1(x1, v), 1e-8)
ineqaulity_constraint!(opt, (x, g) -> const2(x3, v), 1e-8)
equality_constraint!(opt, (x, g)   -> const3(K1, v), 1e-8)

(minf, minx, ret) = optimize(opt, [1.0, 1.0, 1.0])

numevals = opt.numevals # the number of function evaluations

println("got $minf at $minx after $numevals iterations (returned $ret)")