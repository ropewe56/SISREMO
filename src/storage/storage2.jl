include("base.jl")


"""
    compute_storage_level(Load, WWSB, torage_capacity)
"""
function compute_storage_level(Load, WWSB_scaled, storage_capacity)
    nb_steps  = length(WWSB_scaled)
    storage = Storage(storage_capacity, nb_steps)
    storage.fill_level[1] = storage_capacity

    istep = 1
    for istep in 2:nb_steps
        power_step(storage, Load, WWSB_scaled, istep)
    end
    #plt.plot(storage.fill_level)
    #extrema(storage.fill_level)
    storage
end

function get_WWSB_scaled(Load, WWSB, x)
    mean_load = mean(Load)
    mean_wwsb = mean(WWSB)
    scale = (mean_load*x[1]) / mean_wwsb
    WWSB .* scale
end

"""
    compute storage fill level for different combinations of storage_capacity and over_production
"""
function optimize_overproduction(Load, WWSB, storage_capacity)
    
    function ff(x)
        WWSB_scaled_0 = get_WWSB_scaled(Load, WWSB, x)
        storage = compute_storage_level(Load, WWSB_scaled_0, storage_capacity)
        storage_min = minimum(storage.fill_level)
        a = (storage_min - storage_capacity*0.1)^2 + sum(storage.other)
        #@info x[1], storage_min, a 
        a
    end
    #plt.plot(WWSB_scaled_0)

    x = [5.0]
    results = Optim.optimize(ff, x, NelderMead())
    
    op = Optim.minimizer(results)
    WWSB_scaled = get_WWSB_scaled(Load, WWSB, op)
    storage = compute_storage_level(Load, WWSB_scaled, storage_capacity)
    storage_min = minimum(storage.fill_level)
    a = (storage_min - storage_capacity*0.1)^2 + sum(storage.other)
    @info storage_capacity, a, Optim.minimizer(results), Optim.minimum(results)

    storage_capacity, Optim.minimizer(results), Optim.minimum(results)
end

"""
    determine_overproduction(Load, RP)

    determine minimum storage capacity as a function of over production

    Load : detrended Load
    RP : detrended and scaled renewable power
    punit : unit of Load and RP (MW. GW, TW)
"""
function determine_overproduction(Load::Vector{Float64}, WWSB::Vector{Float64}, storage_capacities)
    ops = []
    storage_capacity = storage_capacities[1]
    for storage_capacity in storage_capacities
        push!(ops, optimize_overproduction(Load, WWSB, storage_capacity))
    end
    ops
end

function run_simulation(start_year, end_year)
    par = make_power_parameter(start_year, end_year)
    detrended_power = get_detrended_power(par; load_arrow = true)
    
    storage_capacities = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0] .* uconversion_factor(par.punit, 1u_TW)

    Load = detrended_power.Load
    WWSB = @. detrended_power.Woff + detrended_power.Won + detrended_power.Solar + detrended_power.Bio

    determine_overproduction(Load, WWSB, storage_capacities)

    storage = compute_storage(Load, WWSB, storage_capacities[3], 3.0)
    storage_min = extrema(storage.SF)[1]

    plt.plot(storage.SF)
   
    (storage_min * (1.0 - 0.1))^2 + sum(storage.I7)
end


op = 2.0
st = compute_storage(Load, WWSB, op)

extrema(st.SF)