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