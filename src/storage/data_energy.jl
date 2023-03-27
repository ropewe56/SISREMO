include("data_energy_charts.jl")

using Unitful
import Unitful.uconvert

"""
    convert units
"""
function uconvert(efrom::String, eto::String)
    ut = Dict("MW" => u"MW", "GW" => u"GW", "TW" => u"TW")
    uf = Dict("MW" => 1u"MW", "GW" => 1u"GW", "TW" => 1u"TW")
    ef = uf[efrom]
    et = ut[eto]
    c = Float64(getproperty(uconvert(et, ef), :val))
end

"""
    ise data

    "Load (MW)",                                #   2
    "Residual load (MW)",                       #   3
    "Wind offshore (MW)",                       #   4
    "Wind onshore (MW)",                        #   5
    "Solar (MW)",                               #   6
    "Biomass (MW)",                             #   7
    "Nuclear (MW)",                             #  15

    uts    :: Vector{Int64}    # 1
    Load   :: Vector{Float64}  # 2
    Woff   :: Vector{Float64}  # 3
    Won    :: Vector{Float64}  # 4
    Solar  :: Vector{Float64}  # 5
    Bio    :: Vector{Float64}  # 6
    WoffC  :: Vector{Float64}  # 7
    WonC   :: Vector{Float64}  # 8
    SolarC :: Vector{Float64}  # 9
    Nuclear:: Vector{Float64}  # 10
    dates  :: Vector{DateTime} # 11
"""
function EnergyData_ise(data_dir)
    # load data matrix from hdf5 file
    D = load_ise_as_hdf5(data_dir)
    # Unix time stamps are stored as Float64 in D => to Int
    uts = floor.(Int, D[:,1])
    # convert Unix time stamps to Date objects
    dates = unix2datetime.(uts)
    # return  uts, Load, Wind off, Wind on, Soalr, Biomass, Nuclear, dates
    uts, D[:,2], D[:,4], D[:,5], D[:,6], D[:,7],  D[:,15], dates
end

struct EnergyData
    uts     :: Vector{Int64}    # 1
    Load    :: Vector{Float64}  # 2
    Woff    :: Vector{Float64}  # 3
    Won     :: Vector{Float64}  # 4
    Solar   :: Vector{Float64}  # 5
    Bio     :: Vector{Float64}  # 6
    Nuclear :: Vector{Float64}  # 10
    dates   :: Vector{DateTime} # 11
end


function average_to_hour(X)
    n = size(X,1)
    #m = div(n, 4)
    #Y = zeros(Float64, m)
    #j = 1
    #for i in 1:m
    #    Y[i] = (X[j] + X[j+1] + X[j+2] + X[j+3]) * 0.25
    #    j += 4
    #end
    @. (X[1:4:n-3] + X[2:4:n-2] + X[3:4:n-1] + X[4:4:n]) * 0.25
end

"""
    EnergyData(eunit)

    EnergyData constructor
"""
function EnergyData(data_dir, eunit)
    uts_4, Load_4, Woff_4, Won_4, Solar_4, Bio_4, Nuclear_4, dates_4 = EnergyData_ise(data_dir)
    # energy charts data are in MW, M2 is conversion factor to eunit (MW, GW, TW)

    n = length(Load_4)
    uts     = uts_4[1:4:n]
    dates   = dates_4[1:4:n]
    Load    = average_to_hour(Load_4)
    Woff    = average_to_hour(Woff_4)
    Won     = average_to_hour(Won_4)
    Solar   = average_to_hour(Solar_4)
    Bio     = average_to_hour(Bio_4)
    Nuclear = average_to_hour(Nuclear_4)

    M2 = uconvert("MW", eunit)

    @infoe ("# timesteps =", length(Load), "energy conversion =", M2)

    EnergyData(uts, Load .* M2, Woff .* M2, Won .* M2, Solar .* M2, Bio .* M2, Nuclear .* M2, dates)
end

n = 10
collect(1:4:n)