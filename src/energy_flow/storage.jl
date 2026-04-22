include("energy_flow.jl")
include("powers.jl")
include("costs.jl")

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function minimize_overproduction(powers::Powers, storage_capacity, info)

    function compute_obj(energy_flow, storage_capacity)
        storage_min = minimum(energy_flow.fill_level)
        obj = (storage_min - storage_capacity*0.1)^2 + sum(energy_flow.other)
        obj
    end

    function ff(x)
        Load = powers.Load
        WWSB_scaled = powers.Woff * x[1] + powers.Won * x[1] + powers.Solar * x[1] + powers.Bio
        energy_flow = compute_energy_flow(Load, WWSB_scaled, storage_capacity)
        obj = compute_obj(energy_flow, storage_capacity)
        obj
    end

    x = [5.0]
    results = Optim.optimize(ff, x, NelderMead())
    
    x = Optim.minimizer(results)
    Load = powers.Load
    WWSB_scaled = powers.Woff * x[1] + powers.Won * x[1] + powers.Solar * x[1] + powers.Bio
    energy_flow = compute_energy_flow(Load, WWSB_scaled, storage_capacity)
    
    if info
        storage_min = minimum(energy_flow.fill_level)
        obj = (storage_min - storage_capacity*0.1)^2 + sum(energy_flow.other)
        opmin = Optim.minimizer(results)
        @info "      ", storage_capacity, opmin, obj, Optim.minimum(results)
    end

    Optim.minimizer(results), Optim.minimum(results), energy_flow
end

function minimize_cost(powers::Powers, info)
    info = false
    x = storage_capacity
    function ff(x)
        (op, min_objective, energy_flow) = minimize_overproduction(powers, x[1], info)
        cost1, cost2 = determine_cost(op[1], energy_flow, powers, x[1])
        @info op, x[1], cost1, cost2
        cost1^2 + cost2^2
    end

    x = [2500.0]
    results = Optim.optimize(ff, x, NelderMead())
    Optim.minimizer(results), Optim.minimum(results)
end


"""
    determine_overproduction(Load, RP)

    determine minimum storage capacity as a function of over production

    Load : detrended Load
    RP : detrended and scaled renewable power
    punit : unit of Load and RP (MW. GW, TW)
"""
function determine_overproduction(powers::Powers, storage_capacities, info)
    ops = []
    storage_capacity = storage_capacities[1]
    for storage_capacity in storage_capacities
        push!(ops, minimize_overproduction(powers, storage_capacity, info))
    end
    ops
end

function init_run(start_year, end_year)
    par = make_power_parameter(start_year, end_year)
    detrended_power = get_detrended_power(par; load_arrow = true)

    Load  = detrended_power.Load ;
    Woff  = detrended_power.Woff ;
    Won   = detrended_power.Won  ;
    Solar = detrended_power.Solar;
    Bio   = detrended_power.Bio  ;

    Powers(Load, Woff, Won, Solar, Bio);
end

function run_simulation(start_year=2017, end_year=2024)
    start_year, end_year = 2017, 2024
    info = true
    powers = init_run(start_year, end_year);

    storage_capacities = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0] .* uconversion_factor(par.punit, 1u_TW)
    determine_overproduction(powers, storage_capacities, info)

    energy_flow = compute_storage(powers, storage_capacities[3], 3.0)
    storage_min = extrema(storage.fill_level)[1]

    plt.plot(energy_flow.fill_level)
   
end

