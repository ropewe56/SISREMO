function run()
    bcap = collect(range(10.0, 1.0e3, 3))
    hcap = collect(range(1.0e3, 5.0e4, 3))
    op = 1.4

    cost = create_heatmap(power_data, nhours, p, bcap, hcap, op)


    cost = JLD.load("cost.jld")["cost"]
    extrema(cost)
    cost0 = cost .- minimum(cost)

    plt.imshow(cost0)

    bcap = collect(range(10.0, 1.0e3, 200))
    hcap = collect(range(1.0e3, 5.0e4, 200))

    cost0[1,1]
    X
    X = bcap' .+ 0*hcap
    Y = 0*bcap' .+ hcap

    X = bcap * hcap'
    Y = hcap * bcap'
    fig, ax = plt.subplots(subplot_kw=Dict("projection" =>"3d"))
    surf = ax.plot_surface(X, Y, cost0, linewidth=0, antialiased=false)

    T = Float64
    power_data, ppar = get_power_data();
    #save_detrended_power_data(power_data, "save_detrended_power_data.hdf5")
    power_data = load_detrended_power_data("save_detrended_power_data.hdf5");

    function get_fopt(power_data)
        nhours = length(power_data.dates)
        p = EnergyParameter{T}()
        function fv1(x)
            p.fcall += 1
            compute(x, power_data, nhours, p)[1]
        end
        function fv2(x, p)
            p.fcall += 1
            compute(x, power_data, nhours, p)[1]
        end
        fv1, fv2, p
    end
    fv1, fv2, p = get_fopt(power_data)


    lb = [T(  10.0), T(2.0e3), T(1.3)]
    ub = [T(1000.0), T(5.0e4), T(1.5)]
    u0 = [0.5*(lb[i]+ub[i]) for i in 1:3]

    sol1 = optimize(fv1, lb, ub, u0, NelderMead()) #f_reltol = 0.01))#, x_abstol=1.0e-5, iterations=200))
    u1 = Optim.minimizer(sol1)
    fv2(u1,p)
    u0 = copy(u1)
    u0[3] = 1.3
    function gv2(g, x, p)
        p.gcall += 1
        g[:] = FiniteDifferences.grad(FiniteDifferences.central_fdm(5, 1), fv2, x, p)[1]
    end

    problem = OptimizationProblem(OptimizationFunction(fv2, grad=gv2), u0, p;
                            lb = lb,
                            ub = ub,
                            lcons = nothing,
                            ucons = nothing,
                            sense = nothing)
    sol = solve(problem, Optimization.LBFGS(), maxiters = 100000)#, maxtime = 1000.0)
    #sol = solve(problem, Optimization.NelderMead(), maxiters = 100000)#, maxtime = 1000.0)

    g = zeros(T,3)
    u = copy(sol.u)
    u[3] += 0.01
    fv2(u,p)

    gv2(g, u, p)
end


