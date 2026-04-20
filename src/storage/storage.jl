include("base.jl")

"""
    compute_storage_level(Load, WWSB, torage_capacity)
    Load = power.Load
"""
function compute_storage_level(Load, WWSB_scaled, storage_capacity)
    nb_steps  = length(WWSB_scaled)
    storage = Storage(storage_capacity, nb_steps)
    storage.fill_level[1] = storage_capacity

    istep = 2
    while istep <= nb_steps
        power_step(storage, Load, WWSB_scaled, istep)        
        istep += 1
    end
    
    storage
end

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function minimize_overproduction(power, storage_capacity; info = false)

    function cost(storage, storage_capacity)
        storage_min = minimum(storage.fill_level)
        a = (storage_min - storage_capacity*0.1)^2 + sum(storage.other)
        storage_min, a
    end

    function ff(x)
        Load = power.Load
        #WWSB_scaled_0 = get_WWSB_scaled(power, x[1], 0)
        WWSB_scaled_1 = get_WWSB_scaled(power, x, 1)
        storage = compute_storage_level(Load, WWSB_scaled_1, storage_capacity)
        storage_min, a = cost(storage, storage_capacity)
        #@info x[1], storage_min, a 
        a
    end

    #plt.plot(WWSB_scaled_0)
    #plt.plot(WWSB_scaled2)

    x = [5.0]
    results = Optim.optimize(ff, x, NelderMead())
    
#    op = Optim.minimizer(results)
#    Load = power.Load
#    WWSB_scaled = get_WWSB_scaled(power, op, 1)
#    storage = compute_storage_level(Load, WWSB_scaled, storage_capacity)
#    storage_min = minimum(storage.fill_level)
#    if info
#        a = (storage_min - storage_capacity*0.1)^2 + sum(storage.other)
#        @info "      ", storage_capacity, a, Optim.minimizer(results), Optim.minimum(results)
#    end

    Optim.minimizer(results), Optim.minimum(results)
end

function minimize_cost(power::ScaledPower, storage_capacities; info = true)

    function cost(op, storage_capacity)
        (op + storage_capacity*1.0e-3)
    end

    function ff(x)
        storage_capacity = x[2]
        (op, min_objective) = minimize_overproduction(power, storage_capacity; info = info)
        a = cost(op[1], storage_capacity)
        @info op[1], storage_capacity, min_objective, a
        a
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
function determine_overproduction(power::ScaledPower, storage_capacities; info = true)
    ops = []
    storage_capacity = storage_capacities[1]
    for storage_capacity in storage_capacities
        push!(ops, minimize_overproduction(power, storage_capacity; info = info))
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
    power, storage_capacities = init_run(start_year, end_year);

    determine_overproduction(power, storage_capacities)

    storage = compute_storage(power, storage_capacities[3], 3.0)
    storage_min = extrema(storage.fill_level)[1]

    plt.plot(storage.fill_level)
   
end

start_year, end_year = 2017, 2024
