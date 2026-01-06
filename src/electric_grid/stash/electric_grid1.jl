import PyPlot as plt
plt.pygui(true)

using Optim

include("power_data.jl")

mutable struct Production{T}
    Et   :: Vector{T}
    ΔE   :: Vector{T}
    cMWh :: T
end
function Production(Et::Vector{T}, op::T, cMWh::T=100.0) where T
    Production(copy(Et) .* op, zeros(T,length(Et)), cMWh)
end
function cost_per_GWh(prod::Production{T}, E::T) where T
    prod.cMWh*1.0e3 * E
end

mutable struct Load{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    cost :: T
    E :: T
end
function Load(Et::Vector{T}) where T
    Load(copy(Et), zeros(T, length(Et)), zero(T), zero(T))
end

mutable struct Import{T}
    Et   :: Vector{T}
    ΔE   :: Vector{T}
    cMWh :: T
end
function Import(n, cMWh = 100.0)
    T = typeof(cMWh)
    Import(zeros(T, n), zeros(T, n), cMWh)
end
function cost_per_GWh(imp::Import{T}, E::T) where T
    imp.cMWh*1.0e3 * E
end

mutable struct Export{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    cMWh :: T
end
function Export(n, cMWh=100.0)
    T = typeof(cMWh)
    Export(zeros(T, n), zeros(T, n), cMWh)
end
function cost_per_GWh(ex::Export{T}, E::T) where T
    ex.cMWh*1.0e3 * E
end

abstract type AbstractStorage{T} end

mutable struct Battery{T} <: AbstractStorage{T}
    CAP  :: T
    E    :: Vector{T}
    ΔE   :: Vector{T}
    inP  :: T
    outP :: T
    ηin  :: T
    ηout :: T
    cMWh :: T 
end
function Battery(n::Int64, CAP::T, inP::T, outP::T, ηin::T, ηout::T, cMWh::T = T(100.0); E0::T = zero(T)) where T
    E = zeros(T, n)
    E[1] = E0
    Battery(CAP, E, zeros(T,n), inP, outP, ηin, ηout, cMWh)
end
function cost_per_GWh(bat::Battery{T}, E::T) where T
    bat.cMWh*1.0e3 * E
end

mutable struct Hydrogen{T} <: AbstractStorage{T}
    CAP  :: T
    E    :: Vector{T}
    ΔE   :: Vector{T}
    inP  :: T
    outP :: T
    ηin  :: T
    ηout :: T
    cMWh :: T
end
function Hydrogen(n::Int64, CAP::T, inP::T, outP::T, ηin::T, ηout::T, cMWh::T=200.0; E0 = zero(T)) where T
    E = zeros(T, n)
    E[1] = E0
    Hydrogen(CAP, E, zeros(T,n), inP, outP, ηin, ηout, cMWh)
end
function cost_per_GWh(h2::Hydrogen{T}, E::T) where T
    h2.cMWh*1.0e3 * E
end



function request_source(prod::Production{T}, E, it, Δt) where T
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

function request_source(st::AbstractStorage{T}, E, it, Δt) where T
    function maxenergy_return(e)
        if E/st.ηout <= st.outP*Δt
            st.ΔE[it] = E/st.ηout
        else
            st.ΔE[it] = st.outP*Δt
        end
    end

    if st.E[it] - E >= 0.0
        maxenergy_return(E)
    else
        maxenergy_return(st.E[it])
    end
    if st.E[it] < st.ΔE[it]        
        st.ΔE[it] = st.E[it]
    end
#    if it < 15
#        @info @sprintf("     %6.3f  %6.3f  %6.3f  %6.3f  %s", st.E[it], E/Δt, st.outP, st.ΔE[it], typeof(st))
#    end

    st.ΔE[it] * st.ηout
end
function substract_energy(st::AbstractStorage{T}, it) where T
    st.E[it] -= st.ΔE[it]
end

function request_source(imp::Import{T}, E, it, Δt) where T
    imp.ΔE[it] = E
    E
end
function substract_energy(imp::Import{T}, it) where T
    imp.Et[it] = -imp.ΔE[it]
end

function request_sink(st::AbstractStorage{T}, E, it) where T
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

function request_sink(expp::Export{T}, e, it) where T
    expp.ΔE[it] = e
end
function add_energy(expp::Export{T}, it) where T
    expp.Et[it] += expp.ΔE[it]
end


function consumption(load::Load{T}, sources, it, Δt, io) where T
    Load0 = load.Et[it]
    Load = 0.0
    for source in sources
        Ereq = max(0.0, Load0 - Load)
        Emax = request_source(source, Ereq, it, Δt)
        Load += Emax
        load.cost += cost_per_GWh(source, Emax)
        load.E += Emax
        if (Load0 - Load) <= 0.0
            break
        end
    end
    for source in sources
        substract_energy(source, it)
    end
end

function production(prod::Production{T}, sinks, it, Δt) where T
    Erest = prod.Et[it]
    for sink in sinks
        Emax = request_sink(sink, Erest, it)
        prod.Et[it] -= Emax
        Erest = prod.Et[it]
        if Erest <= 0.0
            break
        end
    end
    for sink in sinks
        add_energy(sink, it)
    end
end

function run_system(load::Load{}, prod::Production{T}, bat::Battery{T}, H2::Hydrogen{T}, 
            impp::Import{T}, expp::Export{T}; Δt=one(T)) where T
    
    io = open("storage_log.txt", "w")
    it = 1
    sources = (prod, bat, H2, impp)
    sinks = (bat, H2, expp)
    for it in eachindex(load.Et)
        if it > 1
            bat.E[it] = bat.E[it-1]
            H2.E[it]  = H2.E[it-1]
        end
        
        #write(io, @sprintf("%5d  %6.3e  %6.3e\n", it, bat.E[it], H2.E[it]))
        
        consumption(load, sources, it, Δt, io)
        production(prod, sinks, it, Δt)
    end
    close(io)
end

T = Float64

power_data, ppar = get_power_data();
nhours = length(power_data.dates)

Base.@kwdef mutable struct EnergyParameter
    prod_cost_factor = 0.7
    prod_CostGWh = 60.0
    bat_CostGWh = 100.0
    H2_CostGWh = 200.0
    import_CostGWh = 80.0
    bat_PowerIn  = 1.0
    bat_PowerOut = 1.0
    H2_PowerIn  = 50.0
    H2_PowerOut = 100.0
    bat_ηin  = 0.9
    bat_ηout = 0.9
    H2_ηin  = 0.7
    H2_ηout = 0.7
    bat_Einit   = 50.0
    H2_Einit   = 1.0e4
end

function compute(x, par)
    bcap, hcap, op = x[1], x[2], x[3]

    load = Load(power_data.Load);

    cMWh = (one(T) + (op-one(T)) * par.prod_cost_factor) * par.prod_CostGWh
    prod = Production(power_data.WWSBPower*op, op, cMWh);
    

    bat  = Battery(nhours,  bcap, par.bat_PowerIn, par.bat_PowerOut, par.bat_ηin, par.bat_ηout, par.bat_CostGWh, E0=par.bat_Einit);
    H2   = Hydrogen(nhours, hcap, par.H2_PowerIn, par.H2_PowerOut, par.H2_ηin, par.H2_ηout, par.H2_CostGWh, E0=par.H2_Einit);

    impp = Import(nhours, par.import_CostGWh)
    expp = Export(nhours)

    run_system(load, prod, bat, H2, impp, expp, Δt=one(T))
    println(x, " ", load.cost)
    abs(load.cost), load, prod, bat, H2, impp, expp
end

par = EnergyParameter()
f(x) = compute(x, par)[1]


initial_x = [ 150.0, 2.0e4, 1.1]
lower     = [   0.0,   0.0, 1.0]
upper     = [1000.0, 5.0e4, 1.5]
sol = optimize(f, lower, upper, initial_x, Fminbox(NelderMead()))#, iterations = 200) #IPNewton())

x = Optim.minimizer(sol)

x = [ 150.0, 2.0e4, 1.0]
f(x)
load.cost, load, prod, bat, H2, impp, expp = compute(x, par);
load.cost
prod.cMWh

plt.plot(power_data.dates, power_data.WWSBPower, label="WWSB")
plt.plot(power_data.dates, power_data.Load, label="L")
plt.legend()

plt.figure()
plt.plot(power_data.dates, bat.E, label="Bat")
plt.legend()
plt.figure()
plt.plot(power_data.dates, H2.E, label="H2")
plt.legend()

plt.figure()
plt.plot(power_data.dates, impp.Et, label="Import")
plt.legend()
plt.figure()
plt.plot(power_data.dates, expp.Et, label="Export")
plt.legend()


plt.figure()
plt.plot(power_data.dates, power_data.Load, label="L")
plt.plot(power_data.dates, prod.ΔE, label="P")
plt.plot(power_data.dates, bat.ΔE, label="B")
plt.plot(power_data.dates, H2.ΔE, label="H")
plt.plot(power_data.dates, impp.ΔE, label="Import")
plt.plot(power_data.dates, expp.ΔE, label="Export")
plt.legend()

500*1.0e9 * 0.2