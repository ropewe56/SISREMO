"""
    Load
"""
mutable struct Load{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    C  :: Vector{T}
    total_cost :: T
    total_energy :: T
end

function Load(Et::Vector{T}) where T
    Load(copy(Et), zeros(T, length(Et)), zeros(T, length(Et)), zero(T), zero(T))
end

function add_cost(load::Load{T}, val) where T
    load.total_cost += val
end

function add_energy(load::Load{T}, E) where T
    load.total_energy += E
end

"""
    Production
"""
mutable struct Production{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    C  :: Vector{T}
    C_GWh :: T
end

function Production(Et::Vector{T}, price_of_MWh::T) where T
    Production(copy(Et), zeros(T,length(Et)), zeros(T,length(Et)), price_of_MWh*1.0e3)
end

function get_energy(prod, it)
    prod.Et[it]
end
function reduce_energy(prod, ΔE, it)
    prod.Et[it] - ΔE
end

function bill_source(prod::Production{T}, ΔE, ΔC, it) where T
    prod.Et[it] -= ΔE
    prod.ΔE[it] += ΔE
    prod.C[it] = ΔC * ΔE
end

function request_source(prod::Production{T}, E::T, it::Int64, Δt::T) where T
    Er = min(E, prod.Et[it])
    Er, prod.C_GWh
end

"""
    Import
"""
mutable struct Import{T}
    Et   :: Vector{T}
    ΔE   :: Vector{T}
    C    :: Vector{T}
    inP  :: T
    C_GWh :: T
end

function Import(n, inP::T, price_of_MWh::T) where T
    Import(zeros(T, n), zeros(T, n), zeros(T, n), inP, price_of_MWh*1.0e3)
end

function bill_source(imp::Import{T}, ΔE, ΔC, it) where T
    imp.Et[it] -= ΔE
    imp.ΔE[it] += ΔE
    imp.C[it] = ΔC * ΔE
end

function request_source(imp::Import{T}, E::T, it, Δt::T) where T
    Er = min(E, imp.inP*Δt)
    Er, imp.C_GWh
end

"""
    Export
"""
mutable struct Export{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    C  :: Vector{T}
    outP :: T
    C_GWh :: T
end

function Export(n, outP, price_of_MWh::T) where T
    Export(zeros(T, n), zeros(T, n), zeros(T, n), outP, price_of_MWh*1.0e3)
end

function bill_sink(expp::Export{T}, ΔE, ΔC, it) where T
    expp.Et[it] -= ΔE
    expp.ΔE[it] += ΔE
    expp.C[it]  -= ΔC * ΔE
end

function request_sink(expp::Export{T}, E::T, it, Δt) where T
    min(E, expp.outP*Δt), expp.C_GWh
end


"""
    Curtailment
"""
mutable struct Curtailment{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    C  :: Vector{T}
    C_GWh :: T
end

function Curtailment(n, price_of_MWh::T) where T
    Curtailment(zeros(T, n), zeros(T, n), zeros(T, n), price_of_MWh*1.0e3)
end

function bill_sink(curt::Curtailment{T}, ΔE, ΔC, it) where T
    curt.Et[it] -= ΔE
    curt.ΔE[it] += ΔE
    curt.C[it] -= ΔC * ΔE
end

function request_sink(curt::Curtailment{T}, E::T, it, Δt) where T
    E, curt.C_GWh
end

"""
    Curtailment
"""
mutable struct ResidualLoad{T}
    Et :: Vector{T}
    ΔE :: Vector{T}
    C  :: Vector{T}
    C_GWh :: T
end

function ResidualLoad(n, price_of_MWh::T) where T
    ResidualLoad(zeros(T, n), zeros(T, n), zeros(T, n), price_of_MWh*1.0e3)
end

function bill_source(res::ResidualLoad{T}, ΔE, ΔC, it) where T
    res.Et[it] -= ΔE
    res.ΔE[it] += ΔE
    res.C[it]  += ΔC * ΔE
end

function request_source(res::ResidualLoad{T}, E::T, it, Δt::T) where T
    E, res.C_GWh
end

function bill_sink(res::ResidualLoad{T}, ΔE, ΔC, it) where T
    res.Et[it] -= ΔE
    res.ΔE[it] += ΔE
    res.C[it]  -= ΔC * ΔE
end

function request_sink(res::ResidualLoad{T}, E::T, it, Δt) where T
    E, res.C_GWh
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
    Ci   :: Vector{T}
    Co   :: Vector{T}
    inP  :: T
    outP :: T
    ηin  :: T
    ηout :: T
    Ci_GWh :: T 
    Co_GWh :: T 
end

function make_battery(n::Int64, CAP, inP, outP, ηin, ηout, Ci_MWh, Co_MWh, E0::T) where T
    E = zeros(T, n)
    E[1] = E0
    Battery(CAP, E, zeros(T,n), zeros(T,n), zeros(T,n), zeros(T,n), inP, outP, ηin, ηout, Ci_MWh*1.0e3, Co_MWh*1.0e3)
end

mutable struct Hydrogen{T} <: AbstractStorage{T}
    CAP  :: T
    E    :: Vector{T}
    ΔEi  :: Vector{T}
    ΔEo  :: Vector{T}
    Ci   :: Vector{T}
    Co   :: Vector{T}
    inP  :: T
    outP :: T
    ηin  :: T
    ηout :: T
    Ci_GWh :: T
    Co_GWh :: T
end

function make_hydrogen(n::Int64, CAP::T, inP::T, outP::T, ηin::T, ηout::T, Ci_MWh::T, Co_MWh::T, E0::T) where T
    E = zeros(T, n)
    E[1] = E0
    Hydrogen(CAP, E, zeros(T,n), zeros(T,n), zeros(T,n), zeros(T,n), inP, outP, ηin, ηout, Ci_MWh*1.0e3, Co_MWh*1.0e3)
end

function bill_source(st::AbstractStorage{T}, ΔE, ΔC, it) where T
    st.E[it] -= ΔE/st.ηout
    st.ΔEo[it] += ΔE
    st.Co[it] = ΔC * ΔE
end

function request_source(st::AbstractStorage{T}, E::T, it, Δt::T) where T
    Eη_max = st.E[it]*st.ηout
    Eo_max = min(Eη_max, st.outP*Δt)
    
    Er = if Eo_max < E
        Eo_max
    else
        E
    end

    Cr = st.Co_GWh

    Er, Cr
end

function bill_sink(st::AbstractStorage{T}, ΔE, ΔC, it) where T
    st.E[it] += ΔE * st.ηin
    st.ΔEi[it] += ΔE
    st.Ci[it] -= ΔC * ΔE
end

function request_sink(st::AbstractStorage{T}, E::T, it, Δt::T) where T
    Eη_max = st.CAP - st.E[it]
    Ei_max = min(Eη_max/st.ηin, st.inP*Δt)

    Er = if Ei_max < E
        Ei_max
    else
        E
    end

    Cr = st.Ci_GWh

    Er, Cr
end


"""
    Load
"""
function consumption(load::Load{T}, sources, it, Δt::T) where T
    ΔL = load.Et[it]

    # sources = (prod, bat, H2, impp, res)
    Er = Vector{T}(undef, length(sources))
    Cr = Vector{T}(undef, length(sources))
    for (i,source) in enumerate(sources)
        Er[i], Cr[i] = request_source(source, load.Et[it], it, Δt)
    end

    ii = sortperm(Cr)

    ΔE = zeros(T, length(sources))
    ΔC = zeros(T, length(sources))
    for i in eachindex(sources)
        j = ii[i]
        if Er[j] > ΔL
            ΔC[j]= Cr[ii[i]]
            ΔE[j] = ΔL
            ΔL -= ΔE[j]
        else
            ΔC[j] = Cr[ii[i]]
            ΔE[j] = Er[j]
            ΔL -= ΔE[j]
        end
        if ΔL <= 0.0
            break
        end
    end
    
    Cmax = maximum(ΔC)
    for i in eachindex(sources)
        j = ii[i]
        load.ΔE[it]+= ΔE[j]
        load.C[it] += ΔE[j] * Cmax
        bill_source(sources[j], ΔE[j], Cmax, it)
    end
end

function production(prod::Production{T}, sinks, it, Δt::T) where T
    ΔP = get_energy(prod, it)

    # sinks = (bat, H2, expp, curt)
    Er = Vector{T}(undef, length(sinks))
    Cr = Vector{T}(undef, length(sinks))
    ΔE = zeros(T, length(sinks))

    for (i,sink) in enumerate(sinks)
        Er[i], Cr[i] = request_sink(sink, ΔP, it, Δt)
    end

    ii = sortperm(Cr, rev=true)

    Cmin = 0.0
    for i in eachindex(sinks)
        j = ii[i]
        if Er[j] > ΔP
            Cmin = Cr[ii[i]]
            ΔE[j] = ΔP
            ΔP = 0.0
        else
            Cmin = Cr[ii[i]]
            ΔE[j] = Er[j]
            ΔP -= Er[j]
        end
        if ΔP <= 0.0
            break
        end
    end

    for i in eachindex(sinks)
        j = ii[i]
        prod.ΔE[it]+= ΔE[j]
        prod.C[it] += ΔE[j] * Cmin
        bill_sink(sinks[j], ΔE[j], Cmin, it)
    end
end

function run_system(load::Load{}, prod::Production{T}, bat::Battery{T}, H2::Hydrogen{T}, 
                    impp::Import{T}, expp::Export{T}, curt::Curtailment{T}, res::ResidualLoad{T}, Δt::T) where T
    
    sources = (prod, bat, H2, impp, res)
    sinks   = (bat, H2, expp, curt)

    it = 1
    for it in eachindex(load.Et)
        if it > 1
            bat.E[it] = bat.E[it-1]
            H2.E[it]  = H2.E[it-1]
        end
        
        consumption(load, sources, it, Δt)
        production(prod, sinks, it, Δt)
    end

end

