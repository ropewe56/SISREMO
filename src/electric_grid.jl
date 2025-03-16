struct NetParameter{T}
    Δt    :: T
    times :: Vector{T}
end

mutable struct Production{T}
    energy_t :: Vector{T}
    ΔE :: T
end

mutable struct Consumption{T}
    energy_t :: Vector{T}
    ΔE :: T
end

mutable struct Import{T}
    energy_t :: Vector{T}
    ΔE :: T
end

mutable struct Export{T}
    energy_t :: Vector{T}
    ΔE :: T
end

abstract type Storage{T} end

mutable struct Battery{T} <: Storage{T}
    capacity :: T
    energy   :: Vector{T}
    inpower  :: T
    outpower :: T
    ηin      :: T
    ηout     :: T
    ΔE :: T
end

mutable struct Hydrogen{T} <: Storage{T}
    capacity :: T
    energy   :: Vector{T}
    inpower  :: T
    outpower :: T
    ηin      :: T
    ηout     :: T
    ΔE :: T
end

function request_energy(prod::Production{T}, e, it, par) where T
    if prod.energy_t[it] >= e
        prod.ΔE = e
    else
        prod.ΔE = prod.energy_t[it]
    end
    prod.ΔE
end

function substract_energy(prod::Production{T}, it) where T
    prod.energy_t[it] -= prod.ΔE
end

function request_energy(st::Storage{T}, e, it, par) where T
    function maxenergy_return(e)
        if e/par.Δt <= st.outpower
            st.ΔE = e
        else
            st.ΔE = st.outpower*par.Δt
        end
    end

    if st.energy[it] - e >= 0.0
        maxenergy_return(e)
    else
        maxenergy_return(st.energy[it])
    end

    st.ΔE * st.ηout
end
function substract_energy(st::Storage{T}, it) where T
    st.energy[it] -= st.ΔE
end

function request_energy(imp::Import{T}, e, it, par) where T
    imp.ΔE = e
    e
end
function substract_energy(imp::Import{T}, it) where T
    imp.energy_t[it] = -imp.ΔE
end

function consumption(load::Consumption{T}, sources, it, par) where T
    le0 = load.energy_t[it]
    le = 0.0
    source = sources[1]
    for source in sources
        ereq = max(0.0, le0 - le)
        emax = request_energy(source, ereq, it, par)
        le += emax
        if (le0 - le) <= 0.0
            break
        end
    end
    for source in sources
        substract_energy(source, it)
    end
end

function production(prod::Production{T}, sinks, it, par) where T
    erest = prod.energy_t[it]
    for sink in sinks
        emax = deliver_energy(sink, erest, it, par)
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

function deliver_energy(st::Storage{T}, e, it, par) where T
    if st.capacity > st.energy[it]
        st.energy[it] += e
        st.ΔE = e * st.ηin
    else
        st.ΔE = (st.capacity - st.energy[it]) * st.ηin
    end
end
function add_energy(st::Storage{T}, it) where T
    st.energy[it] = st.ΔE
end

function deliver_energy(expp::Export{T}, e, it, par) where T
    expp.ΔE = e
end
function add_energy(expp::Export{T}, it) where T
    expp.energy_t[it] += expp.ΔE
end

function run_system(load::Consumption{}, prod::Production{T}, bat::Storage{T}, h2::Storage{T}, 
            impp::Import{T}, expp::Export{T}, par) where T
    
    it = 1
    sources = (prod, bat, h2, impp)
    sinks = (bat, h2, expp)
    for it in eachindex(par.times)
        if it > 1
            bat.energy[it] = bat.energy[it-1]
            h2.energy[it] = h2.energy[it-1]
        end
        consumption(load, sources, it, par)
        production(prod, sinks, it, par)
    end
end


function setup_system(par, rl, rp, ::Type{T}=Float64) where T
    
    le = (rl .+ T(4.0))
    pe = (rp .+ T(0.1))
    
    pen = pe * sum(le) / sum(pe)
    pe0 = copy(pen)
    
    load = Consumption(le, 0.0)
    prod = Production(pen, 0.0)

    etot = sum(pe0)
    bat  = Battery(etot*0.1, zeros(T, nhours), 1.0, 1.0, 0.8, 0.8, 0.0)
    h2   = Hydrogen(etot*0.5, zeros(T, nhours), 1.0, 1.0, 0.6, 0.6, 0.0)

    impp = Import(zeros(T, nhours), 0.0)
    expp = Export(zeros(T, nhours), 0.0)

    run_system(load, prod, bat, h2, impp, expp, par)

    prod.energy_t = pe0
    load, prod, bat, h2, impp, expp
end
T = Float64
nhours = 24*365
par = NetParameter(1.0, ones(T, nhours))
rl = rand(T, nhours)
rp = rand(T, nhours)

load, prod, bat, h2, impp, expp = setup_system(par, rl, rp)
using Printf
for i in 1:10
    @printf("%5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n", load.energy_t[i], prod.energy_t[i], bat.energy[i], h2.energy[i], impp.energy_t[i], expp.energy_t[i])
end