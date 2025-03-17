import PyPlot as plt
plt.pygui(true)

include("power_data.jl")

mutable struct Production{T}
    energy_t :: Vector{T}
    ΔE       :: Vector{T}
end
function Production(energy_t::Vector{T}) where T
    Production(copy(energy_t), zeros(T, length(energy_t)))
end

mutable struct Consumption{T}
    energy_t :: Vector{T}
    ΔE       :: Vector{T}
end
function Consumption(energy_t::Vector{T}) where T
    Consumption(copy(energy_t), zeros(T, length(energy_t)))
end

mutable struct Import{T}
    energy_t :: Vector{T}
    ΔE       :: Vector{T}
end
function Import(n, ::Type{T}) where T
    Import(zeros(T, n), zeros(T, n))
end

mutable struct Export{T}
    energy_t :: Vector{T}
    ΔE :: Vector{T}
end
function Export(n, ::Type{T}) where T
    Export(zeros(T, n), zeros(T, n))
end

abstract type AbstractStorage{T} end

mutable struct Battery{T} <: AbstractStorage{T}
    capacity :: T
    energy   :: Vector{T}
    ΔE       :: Vector{T}
    inpower  :: T
    outpower :: T
    ηin      :: T
    ηout     :: T
end
function Battery(n::Int64, c::T, ip::T, op::T, ηin::T, ηout::T; E0::T = zero(T)) where T
    e = zeros(T, n)
    e[1] = E0
    Battery(c, e, zeros(T,n), ip, op, ηin, ηout)
end

mutable struct Hydrogen{T} <: AbstractStorage{T}
    capacity :: T
    energy   :: Vector{T}
    ΔE       :: Vector{T}
    inpower  :: T
    outpower :: T
    ηin      :: T
    ηout     :: T
end
function Hydrogen(n::Int64, c::T, ip::T, op::T, ηin::T, ηout::T; E0 = zero(T)) where T
    e = zeros(T, n)
    e[1] = E0
    Hydrogen(c, e, zeros(T,n), ip, op, ηin, ηout)
end



function request_energy(prod::Production{T}, e, it, Δt) where T
    if prod.energy_t[it] >= e
        prod.ΔE[it] = e
    else
        prod.ΔE[it] = prod.energy_t[it]
    end
    prod.ΔE[it]
end
function substract_energy(prod::Production{T}, it) where T
    prod.energy_t[it] -= prod.ΔE[it]
end

function request_energy(st::AbstractStorage{T}, e, it, Δt) where T
    function maxenergy_return(e)
        if e/Δt <= st.outpower
            st.ΔE[it] = e
        else
            st.ΔE[it] = st.outpower*Δt
        end
    end

    if st.energy[it] - e >= 0.0
        maxenergy_return(e)
    else
        maxenergy_return(st.energy[it])
    end

    st.ΔE[it] * st.ηout
end
function substract_energy(st::AbstractStorage{T}, it) where T
    st.energy[it] -= st.ΔE[it]
end

function request_energy(imp::Import{T}, e, it, Δt) where T
    imp.ΔE[it] = e
    e
end
function substract_energy(imp::Import{T}, it) where T
    imp.energy_t[it] = -imp.ΔE[it]
end

function deliver_energy(st::AbstractStorage{T}, e, it) where T
    if st.capacity > st.energy[it]
        st.energy[it] += e
        st.ΔE[it] = e * st.ηin
    else
        st.ΔE[it] = (st.capacity - st.energy[it]) * st.ηin
    end
end
function add_energy(st::AbstractStorage{T}, it) where T
    st.energy[it] = st.ΔE[it]
end

function deliver_energy(expp::Export{T}, e, it) where T
    expp.ΔE[it] = e
end
function add_energy(expp::Export{T}, it) where T
    expp.energy_t[it] += expp.ΔE[it]
end


function consumption(load::Consumption{T}, sources, it, Δt) where T
    le0 = load.energy_t[it]
    le = 0.0
    for source in sources
        ereq = max(0.0, le0 - le)
        emax = request_energy(source, ereq, it, Δt)
        le += emax
        if (le0 - le) <= 0.0
            break
        end
    end
    for source in sources
        substract_energy(source, it)
    end
end

function production(prod::Production{T}, sinks, it, Δt) where T
    erest = prod.energy_t[it]
    for sink in sinks
        emax = deliver_energy(sink, erest, it)
        prod.energy_t[it] -= emax
        erest = prod.energy_t[it]
        if erest <= 0.0
            break
        end
    end
    for sink in sinks
        add_energy(sink, it)
    end
end

function run_system(load::Consumption{}, prod::Production{T}, bat::Battery{T}, h2::Hydrogen{T}, 
            impp::Import{T}, expp::Export{T}; Δt=one(T)) where T
    
    it = 1
    sources = (prod, bat, h2, impp)
    sinks = (bat, h2, expp)
    for it in eachindex(load.energy_t)
        if it > 1
            bat.energy[it] = bat.energy[it-1]
            h2.energy[it] = h2.energy[it-1]
        end
        consumption(load, sources, it, Δt)
        production(prod, sinks, it, Δt)
    end
end

T = Float64

power_data, ppar = get_power_data();
nhours = length(power_data.dates)

load = Consumption(power_data.Load);
prod = Production(power_data.WWSBPower);

bat  = Battery(nhours,  150.0,  1.0, 1.0, 0.8, 0.8, E0=0.0);
h2   = Hydrogen(nhours, 50.0e3, 10.0, 10.0, 0.6, 0.6, E0=0.0);

impp = Import(nhours, T)
expp = Export(nhours, T)

run_system(load, prod, bat, h2, impp, expp, Δt=one(T))

plt.plot(power_data.dates, power_data.WWSBPower, label="WWSB")
plt.plot(power_data.dates, power_data.Load, label="L")
plt.legend()

plt.figure()
plt.plot(power_data.dates, bat.energy, label="Bat")
plt.plot(power_data.dates, h2.energy, label="H2")
plt.legend()

plt.figure()
plt.plot(power_data.dates, impp.energy_t, label="Import")
plt.plot(power_data.dates, expp.energy_t, label="Export")
plt.legend()


plt.figure()
plt.plot(power_data.dates, power_data.Load, label="L")
plt.plot(power_data.dates, prod.ΔE, label="P")
plt.plot(power_data.dates, bat.ΔE, label="B")
plt.plot(power_data.dates, h2.ΔE, label="H")
plt.plot(power_data.dates, impp.ΔE, label="Import")
plt.plot(power_data.dates, expp.ΔE, label="Export")
plt.legend()

