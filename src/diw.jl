using CSV
using DataFrames
using RWLogger
using ForwardDiff
using NLopt

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
    L :: Vector{Float64}
    M :: Vector{Float64}
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
    L = Vector{Float64}(undef, 0)
    M = Vector{Float64}(undef, 0)
    Vectors(B, C, D, E, F, G, H, I, J, K, L, M)
end

function load_data()
    ii = Dict(2012 => 3, 2013 => 5, 2014 => 7, 2015 => 9, 2016 => 11)
    i = ii[2016]
    v = Vectors()

    root = "/home/wester/Projects/Julia/Private/Energy.jl/jupyter/Energy/DIW/eer-d-18-00186_code_data/model_and_data_spreadsheet_tool/"
    csv = "spreadsheet_tool_1.csv"
    data = CSV.read(joinpath(root, csv), DataFrame, header=4, delim=';', decimal=',');
    v.B = data[:,i];
    v.C = data[:,i+1];
    v
end

function power_flow(x, v)

    EfficiencyLoading     = 0.81
    EfficiencyDischarging = 0.926

    C2 = x[1] # Renewable capacity  [GW]
    C3 = x[2] # Curtailment threshold [GW]
    K1 = x[3] # storage energy level v.K[1] [GWh]

    # v.B demand
    # v.C renewable availbility factor

    v.D = v.C .* C2                                 # renewable availability
    v.F = min.(v.D, v.B)                            # renewable generation
    v.E = v.D - v.F                                 # renewable surplus = renewable availability - renewable generation
    v.G = v.B - v.F                                 # residual demand   = demand - renewable generation
    v.H = min.(C3 .- v.E, 0.0)                      # curtailment       = Curtailment threshold - renewable surplus
    v.I = v.E - v.H;                                # storage loading   = min(renewable surplus - curtailment, 0.0)

    n = length(v.D)
    v.J = zeros(Float64, n)                         # storage discharging = min(residual demand, storage energy level * eta_out)
    v.J[1] = 0.0

    v.K = zeros(Float64, n)                         # storage energy level += storage loading * eta_in -  storage discharging / eta_out
    v.K[1] = K1

    for i in 2:n
        v.J[i] = min(v.G[i], v.K[i-1] .* EfficiencyDischarging);
        v.K[i] = v.K[i-1] + v.I[i] .* EfficiencyLoading - v.J[i] ./ EfficiencyDischarging;
    end

    v.L = map(x -> max(x, 0.0), v.G .- v.J)          # conventional generation
    v.M = v.B .- v.D;                                # residual demand (for RLDCs)

    maximum(v.K)
end

function constraint1(x1, v)
    # x1 = Storage energy in GWh
    # (renewable generation + storage discharging) / renewable availability
    c1 = (sum(v.F) + sum(v.J)) / sum(v.B) - x1 # = 0
    c1
end

function constraint2(x2, v)
    # x2 = Renewable curtailment in per mille or percent
    # curtailment / renewable availability
    c2 = sum(v.H) / sum(v.D) - x2              # <= 0
    c2
end

function constraint3(K1, v)
    c3 = K1 - v.K[end]                         # = 0
    c3
end

function cconstraint1(x, g, x1, v)
    function f(x)
        r = power_flow(x, v)
        (sum(v.F) + sum(v.J)) / sum(v.B) - x1
    end
    if length(g) > 0
        ForwardDiff.gradient!(g, f, x)
    end
    f(x)
end

function cconstraint2(x, g, x2, v)
    function f(x)
        r = power_flow(x, v)
        sum(v.H) / sum(v.D) - x2
    end
    if length(g) > 0
        ForwardDiff.gradient!(g, f, x)
    end
    f(x)
end

function cconstraint3(x, g, K1, v)
    function f(x)
        r = power_flow(x, v)
        K1 - v.K[end]
    end
    if length(g) > 0
        ForwardDiff.gradient!(g, f, x)
    end
    f(x)
end

function myfunc(x, g, v)
    f(x) = power_flow(x, v)
    if length(g) > 0
        ForwardDiff.gradient!(g, f, x)
    end
    return power_flow(x, v)
end

function optimize_power_flow(x0, x1, x2, v)

    # LN_COBYLA, https://github.com/JuliaOpt/NLopt.jl/blob/master/src/NLopt.jl
    opt = Opt(:LN_COBYLA, 3)
    opt.lower_bounds = [0.0, 0.0, 0.0]
    opt.xtol_rel = 1.0e-8

    min_objective!(opt, (x,g) -> myfunc(x,g,v))

    equality_constraint!(opt,   (x, g) -> cconstraint1(x, g, x1,   v), 1.0e-8)
    inequality_constraint!(opt, (x, g) -> cconstraint2(x, g, x2,   v), 1.0e-8)
    equality_constraint!(opt,   (x, g) -> cconstraint3(x, g, x[3], v), 1.0e-8)

    (minf, minx, ret) = optimize(opt, x0)

    numevals = opt.numevals # the number of function evaluations

    minf, minx, ret, numevals
end

# Renewable share
x1 = collect(LinRange(2.0e-2, 1.0, 10))
# Renewable curtailment
x2 = collect(LinRange(0.0, 0.5, 10))

v = load_data();

# C2 = x[1] # Renewable capacity [GW]
# C3 = x[2] # Curtailment threshold [GW]
# K1 = x[3] # storage energy level v.K[1]
x0 = [200.0, 20.0, 1140.0]

minf, minx, ret, numevals = optimize_power_flow(x0, x1[end], x2[end], v)

constraint1(x1[end], v)
constraint1(x2[end], v)
constraint3(v.K[1],  v)

#v.K
sum(v.L)
sum(v.I)
sum(v.H)

xx1 = (sum(v.F) + sum(v.J)) / sum(v.B);
A = (v.F + v.J) / v.B;
minimum(A), maximum(A)