include("power_data.jl")

"""
    Load
"""
mutable struct Load{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    total_cost :: T
    total_energy :: T
end

function Load(Et::Vector{T}) where T
    Load(copy(Et), zeros(eltype(Et), length(Et)), zero(T), zero(T))
end

function add_cost(load::Load{T}, val) where T
    load.total_cost += val
end

function store_in(load::Load{T}, E) where T
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

function request_source(prod::Production{T}, E::T, it::Int64, Δt::T) where T
    if prod.Et[it] >= E
        prod.ΔE[it] = E
    else
        prod.ΔE[it] = prod.Et[it]
    end
    prod.ΔE[it]
end
function storage_out(prod::Production{T}, it) where T
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

function request_source(imp::Import{T}, E::T, it, Δt::T) where T
    imp.ΔE[it] = E
    E
end

function storage_out(imp::Import{T}, it) where T
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

function request_sink(expp::Export{T}, e::T, it, Δt::T) where T
    expp.ΔE[it] = e
end

function store_in(expp::Export{T}, it) where T
    expp.Et[it] += expp.ΔE[it]
end

"""
    Storage
"""
abstract type AbstractStorage{T} end

mutable struct Battery{T} <: AbstractStorage{T}
    CAP  :: T
    E    :: Vector{T}
    ΔEi  :: Vector{T}
    ΔEo  :: Vector{T}
    inP  :: T
    outP :: T
    ηin  :: T
    ηout :: T
    price_of_MWh :: T 
end

function make_battery(n::Int64, CAP, inP, outP, ηin, ηout, price_of_MWh, E0::T) where T
    E = zeros(T, n)
    E[1] = E0
    Battery(CAP, E, zeros(T,n), zeros(T,n), inP, outP, ηin, ηout, price_of_MWh)
end

function energy_cost(bat::Battery{T}, E::T) where T
    bat.price_of_MWh*1.0e3 * E
end

mutable struct Hydrogen{T} <: AbstractStorage{T}
    CAP  :: T
    E    :: Vector{T}
    ΔEi  :: Vector{T}
    ΔEo  :: Vector{T}
    inP  :: T
    outP :: T
    ηin  :: T
    ηout :: T
    price_of_MWh :: T
end

function make_hydrogen(n::Int64, CAP::T, inP::T, outP::T, ηin::T, ηout::T, price_of_MWh::T, E0::T) where T
    E = zeros(T, n)
    E[1] = E0
    Hydrogen(CAP, E, zeros(T,n), zeros(T,n), inP, outP, ηin, ηout, price_of_MWh)
end

function energy_cost(h2::Hydrogen{T}, E::T) where T
    h2.price_of_MWh*1.0e3 * E
end

function request_source(st::AbstractStorage{T}, E::T, it, Δt::T) where T
    function maxenergy_out(E)
        if E <= st.outP*Δt
            st.ΔEo[it] = E/st.ηout
        else
            st.ΔEo[it] = st.outP*Δt/st.ηout
        end
    end

    if st.E[it] - E >= zero(T)
        maxenergy_out(E)
    else
        maxenergy_out(st.E[it])
    end

    if st.E[it] < st.ΔEo[it]        
        st.ΔEo[it] = st.E[it]
    end

    st.ΔEo[it] * st.ηout
end

function storage_out(st::AbstractStorage{T}, it) where T
    st.E[it] -= st.ΔEo[it]
end

function request_sink(st::AbstractStorage{T}, E::T, it, Δt::T) where T
    function maxenergy_in(E)
        if E <= st.inP*Δt
            st.ΔEi[it] = E
        else
            st.ΔEi[it] = st.inP*Δt
        end
    end

    ΔEmax_in = st.CAP - st.E[it]
    if ΔEmax_in > E
        maxenergy_in(E)
    else
        maxenergy_in(ΔEmax_in)
    end

    st.ΔEi[it]
end

function store_in(st::AbstractStorage{T}, it) where T
    st.ΔEi[it] = st.ΔEi[it] * st.ηin
    st.E[it] += st.ΔEi[it]
end

"""
    Load
"""
function consumption(load::Load{T}, sources, it, Δt::T) where T
    load0 = load.Et[it]
    load = 0.0
    # sources = (prod, bat, H2, impp)
    for source in sources
        Edemand = max(0.0, load0 - load)
        Emax = request_source(source, Edemand, it, Δt)
        Load += Emax
        
        add_cost(load, energy_cost(source, Emax))
        store_in(load, Emax)

        if (Load0 - Load) <= 0.0
            break
        end
    end

    for source in sources
        storage_out(source, it)
    end
end

function production(prod::Production{T}, sinks, it, Δt::T) where T
    E_remaining = get_energy(prod, it)

    for sink in sinks
        Emax = request_sink(sink, E_remaining, it, Δt)

        substract_from_remaining_energy(prod, it, Emax)
        E_remaining = get_energy(prod, it)

        if E_remaining <= zero(T)
            break
        end
    end

    for sink in sinks
        store_in(sink, it)
    end
end

function run_system(load::Load{}, prod::Production{T}, bat::Battery{T}, H2::Hydrogen{T}, 
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

