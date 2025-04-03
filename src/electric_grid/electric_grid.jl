include("power_data.jl")

"""
    Consumption
"""
mutable struct Consumption{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    total_cost :: T
    total_energy :: T
end

function Consumption(Et::Vector{T}) where T
    Consumption(copy(Et), zeros(T, length(Et)), zero(T), zero(T))
end

function add_cost(load::Consumption{T}, val) where T
    load.total_cost += val
end

function add_energy(load::Consumption{T}, E) where T
    load.total_energy += E
end

"""
    Production
"""
mutable struct Production{T}
    Et   :: Vector{T}
    ΔE   :: Vector{T}
    price_of_MWh :: T
end

function Production(Et::Vector{T}, op::T, price_of_MWh::T) where T
    Production(copy(Et) .* op, zeros(T,length(Et)), price_of_MWh)
end

function energy_cost(prod::Production{T}, E::T) where T
    prod.price_of_MWh*1.0e3 * E
end

function request_energy(prod::Production{T}, E::T, it::Int64, Δt::T) where T
    if prod.Et[it] >= E
        prod.ΔE[it] = E
    else
        prod.ΔE[it] = prod.Et[it]
    end
    prod.ΔE[it]
end
function substract_energy(prod::Production{T}, it) where T
    prod.Et[it] -= prod.ΔE[it]
end

function substract_from_remaining_energy(prod::Production{T}, it, E::T) where T
    prod.Et[it] -= E
end
function get_energy(prod::Production{T}, it) where T
    prod.Et[it]
end


"""
    Import
"""
mutable struct Import{T}
    Et   :: Vector{T}
    ΔE   :: Vector{T}
    price_of_MWh :: T
end

function Import(n, price_of_MWh)
    T = typeof(price_of_MWh)
    Import(zeros(T, n), zeros(T, n), price_of_MWh)
end

function energy_cost(imp::Import{T}, E::T) where T
    imp.price_of_MWh*1.0e3 * E
end

function request_energy(imp::Import{T}, E::T, it, Δt::T) where T
    imp.ΔE[it] = E
    E
end

function substract_energy(imp::Import{T}, it) where T
    imp.Et[it] = -imp.ΔE[it]
end

"""
    Export
"""
mutable struct Export{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    price_of_MWh :: T
end

function Export(n, price_of_MWh::T) where T
    Export(zeros(T, n), zeros(T, n), price_of_MWh)
end

function energy_cost(ex::Export{T}, E::T) where T
    ex.price_of_MWh*1.0e3 * E
end

function energy_to_sink(expp::Export{T}, e::T, it) where T
    expp.ΔE[it] = e
end

function add_energy(expp::Export{T}, it) where T
    expp.Et[it] += expp.ΔE[it]
end

"""
    Storage
"""
abstract type AbstractStorage{T} end

mutable struct Battery{T} <: AbstractStorage{T}
    CAP  :: T
    E    :: Vector{T}
    ΔE   :: Vector{T}
    inP  :: T
    outP :: T
    ηin  :: T
    ηout :: T
    price_of_MWh :: T 
end

function make_battery(n::Int64, CAP, inP, outP, ηin, ηout, price_of_MWh, E0::T) where T
    E = zeros(T, n)
    E[1] = E0
    Battery(CAP, E, zeros(T,n), inP, outP, ηin, ηout, price_of_MWh)
end

function energy_cost(bat::Battery{T}, E::T) where T
    bat.price_of_MWh*1.0e3 * E
end

mutable struct Hydrogen{T} <: AbstractStorage{T}
    CAP  :: T
    E    :: Vector{T}
    ΔE   :: Vector{T}
    inP  :: T
    outP :: T
    ηin  :: T
    ηout :: T
    price_of_MWh :: T
end

function make_hydrogen(n::Int64, CAP::T, inP::T, outP::T, ηin::T, ηout::T, price_of_MWh::T, E0::T) where T
    E = zeros(T, n)
    E[1] = E0
    Hydrogen(CAP, E, zeros(T,n), inP, outP, ηin, ηout, price_of_MWh)
end

function energy_cost(h2::Hydrogen{T}, E::T) where T
    h2.price_of_MWh*1.0e3 * E
end

function request_energy(st::AbstractStorage{T}, E::T, it, Δt::T) where T
    function maxenergy_return(e)
        if E/st.ηout <= st.outP*Δt
            st.ΔE[it] = E/st.ηout
        else
            st.ΔE[it] = st.outP*Δt
        end
    end

    if st.E[it] - E >= zero(T)
        maxenergy_return(E)
    else
        maxenergy_return(st.E[it])
    end
    if st.E[it] < st.ΔE[it]        
        st.ΔE[it] = st.E[it]
    end

    st.ΔE[it] * st.ηout
end

function substract_energy(st::AbstractStorage{T}, it) where T
    st.E[it] -= st.ΔE[it]
end

function energy_to_sink(st::AbstractStorage{T}, E::T, it) where T
    ΔEmax = st.CAP - st.E[it]
    if ΔEmax > E
        st.ΔE[it] = E
    else
        st.ΔE[it] = ΔEmax
    end
    st.ΔE[it]
end
function add_energy(st::AbstractStorage{T}, it) where T
    st.E[it] += st.ΔE[it] * st.ηin
end

"""
    Consumption
"""
function consumption(load::Consumption{T}, sources, it, Δt::T) where T
    Load0 = load.Et[it]
    Load = 0.0
    # sources = (prod, bat, H2, impp)
    for source in sources
        Edemand = max(0.0, Load0 - Load)
        Emax = request_energy(source, Edemand, it, Δt)
        Load += Emax
        
        add_cost(load, energy_cost(source, Emax))
        add_energy(load, Emax)

        if (Load0 - Load) <= 0.0
            break
        end
    end

    for source in sources
        substract_energy(source, it)
    end
end

function production(prod::Production{T}, sinks, it, Δt::T) where T
    E_remaining = get_energy(prod, it)

    for sink in sinks
        Emax = energy_to_sink(sink, E_remaining, it)

        substract_from_remaining_energy(prod, it, Emax)
        E_remaining = get_energy(prod, it)

        if E_remaining <= zero(T)
            break
        end
    end

    for sink in sinks
        add_energy(sink, it)
    end
end

function run_system(load::Consumption{}, prod::Production{T}, bat::Battery{T}, H2::Hydrogen{T}, 
                    impp::Import{T}, expp::Export{T}, Δt::T) where T
    
    sources = (prod, bat, H2, impp)
    sinks   = (bat, H2, expp)

    for it in eachindex(load.Et)
        if it > 1
            bat.E[it] = bat.E[it-1]
            H2.E[it]  = H2.E[it-1]
        end
        
        consumption(load, sources, it, Δt)
        production(prod, sinks, it, Δt)
    end

end

Base.@kwdef mutable struct EnergyParameter{T}
    pcfactor :: T = T(0.7)
    pcGWh    :: T = T(60.0)
    bcGWh    :: T = T(100.0)
    hcGWh    :: T = T(200.0)
    icGWh    :: T = T(150.0)
    ecGWh    :: T = T(50.0)
    bPin     :: T = T(1.0)
    bPout    :: T = T(1.0)
    hPin     :: T = T(50.0)
    hPout    :: T = T(100.0)
    bηin     :: T = T(0.9)
    bηout    :: T = T(0.9)
    hηin     :: T = T(0.7)
    hηout    :: T = T(0.7)
    bE0      :: T = T(50.0)
    hE0      :: T = T(1.0e4)
end

