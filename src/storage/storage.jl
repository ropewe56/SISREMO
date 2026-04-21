include("energy_flow.jl")

"""
    compute_energy_flow(Load, WWSB, torage_capacity)
    Load = power.Load
"""
function compute_energy_flow(Load, WWSB_scaled, storage_capacity)
    nb_steps  = length(WWSB_scaled)
    energy_flow = EnergyFlow(storage_capacity, nb_steps)
    energy_flow.fill_level[1] = storage_capacity

    istep = 2
    while istep <= nb_steps
        power_step(energy_flow, Load, WWSB_scaled, istep)        
        istep += 1
    end
    
    energy_flow
end

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function minimize_overproduction(power, storage_capacity, info)

    function compute_obj(energy_flow, storage_capacity)
        storage_min = minimum(energy_flow.fill_level)
        obj = (storage_min - storage_capacity*0.1)^2 + sum(energy_flow.other)
        obj
    end

    function ff(x)
        Load = power.Load
        #WWSB_scaled_0 = get_WWSB_scaled(power, x[1], 0)
        WWSB_scaled_1 = get_WWSB_scaled(power, x, 1)
        energy_flow = compute_energy_flow(Load, WWSB_scaled_1, storage_capacity)
        obj = compute_obj(energy_flow, storage_capacity)
        obj
    end

    x = [5.0]
    results = Optim.optimize(ff, x, NelderMead())
    
    op = Optim.minimizer(results)
    Load = power.Load
    WWSB_scaled = get_WWSB_scaled(power, op, 1)
    energy_flow = compute_energy_flow(Load, WWSB_scaled, storage_capacity)
    
    if info
        storage_min = minimum(energy_flow.fill_level)
        obj = (storage_min - storage_capacity*0.1)^2 + sum(energy_flow.other)
        opmin = Optim.minimizer(results)
        @info "      ", storage_capacity, opmin, obj, Optim.minimum(results)
    end

    Optim.minimizer(results), Optim.minimum(results), energy_flow
end

function minimize_cost(power::ScaledPower, storage_capacities, info)

    function determine_cost(op, storage_capacity, energy_flow)
        (op + storage_capacity*1.0e-3)
    end

    function ff(x)
        storage_capacity = x[2]
        (op, min_objective, energy_flow) = minimize_overproduction(power, storage_capacity, info)
        cost = determine_cost(op[1], storage_capacity, energy_flow)
        @info op[1], storage_capacity, min_objective, a
        cost
    end

    x = [5.0, storage_capacities[3]]
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
function determine_overproduction(power::ScaledPower, storage_capacities, info)
    ops = []
    storage_capacity = storage_capacities[1]
    for storage_capacity in storage_capacities
        push!(ops, minimize_overproduction(power, storage_capacity, info))
    end
    ops
end

function init_run(start_year, end_year)
    par = make_power_parameter(start_year, end_year)
    detrended_power = get_detrended_power(par; load_arrow = true)
    detrended_power.WWSBPower

    Load  = detrended_power.Load ;
    Woff  = detrended_power.Woff ;
    Won   = detrended_power.Won  ;
    Solar = detrended_power.Solar;
    Bio   = detrended_power.Bio  ;
    @. Woff + Won + Solar + Bio

    storage_capacities = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0] .* uconversion_factor(par.punit, 1u_TW)
    power = ScaledPower(Load, Woff, Won, Solar, Bio);
    power.WWSB
    power, storage_capacities
end

function run_simulation(start_year=2017, end_year=2024)
    start_year, end_year = 2017, 2024
    info = true
    power, storage_capacities = init_run(start_year, end_year);

    determine_overproduction(power, storage_capacities, info)

    energy_flow = compute_storage(power, storage_capacities[3], 3.0)
    storage_min = extrema(storage.fill_level)[1]

    plt.plot(energy_flow.fill_level)
   
end

