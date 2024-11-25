include("ise_energy_charts.jl")

using Dates
using PhysConst.UnitConst

function uconversion_factor(eto, efrom)
    @infoe eto, efrom
    Float64(getproperty(uconvert(eto, efrom), :val))
end

function average_to_hour(X)
    n = length(X)
    @. (X[1:4:n-3] + X[2:4:n-2] + X[3:4:n-1] + X[4:4:n]) * 0.25
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
function load_ise_public_power(hdf5_path)
    # load data matrix from hdf5 file
    D, prodtypes = load_ise_data_from_hdf5(hdf5_path)
    id = Dict()
    for (i,pr) in enumerate(prodtypes)
        id[pr] = i+1
    end
    # Unix time stamps are stored as Float64 in D => to Int
    uts = floor.(Int, D[:,1])
    # convert Unix time stamps to Date objects
    dates = unix2datetime.(uts)
    i1 = id["Load"]
    i2 = id["Wind offshore"]
    i3 = id["Wind onshore"]
    i4 = id["Solar"]
    i5 = id["Biomass"]
    i6 = id["Nuclear"]
    # return  uts, 
    uts, D[:,i1], D[:,i2], D[:,i3], D[:,i4], D[:,i5],  D[:,i6], dates
end

struct PowerData
    dates   :: Vector{DateTime} # 1
    uts     :: Vector{Int64}    # 2
    Load    :: Vector{Float64}  # 3
    Woff    :: Vector{Float64}  # 4
    Won     :: Vector{Float64}  # 5
    Solar   :: Vector{Float64}  # 6
    Bio     :: Vector{Float64}  # 7
    Nuclear :: Vector{Float64}  # 8
end

"""
    EnergyData(eunit)

    EnergyData constructor
"""
function PowerData(data_dir, punit, start_year, end_year, scale_Bio=1.0)

    hdf5_path = joinpath(data_dir, @sprintf("public_power_%d-%d.hdf5", start_year, end_year))
    uts_4, Load_4, Woff_4, Won_4, Solar_4, Bio_4, Nuclear_4, dates_4 = load_ise_public_power(hdf5_path)
    
    # energy charts data are in MW, MW_to_unit is conversion factor to eunit (MW, GW, TW)

    MW_to_unit = uconversion_factor(u_MW, 1.0*punit)
    Load    = average_to_hour(Load_4)     .* MW_to_unit
    Woff    = average_to_hour(Woff_4)     .* MW_to_unit
    Won     = average_to_hour(Won_4)      .* MW_to_unit
    Solar   = average_to_hour(Solar_4)    .* MW_to_unit
    Bio     = average_to_hour(Bio_4)      .* scale_Bio .* MW_to_unit
    Nuclear = average_to_hour(Nuclear_4)  .* MW_to_unit

    nl  = length(Load)
    uts = uts_4[1:4:end][1:nl]
    dates = Dates.unix2datetime.(uts)

    @infoe length(uts), length(dates), length(Load)

    @infoe @sprintf("# timesteps = %d, length(Load) = %d, energy conversion = %d", length(uts), length(Load), MW_to_unit)

    PowerData(dates, uts, Load, Woff, Won, Solar, Bio, Nuclear)
end

struct InstalledPower
    dates   :: Vector{DateTime} # 1
    uts     :: Vector{Int64}    # 2
    Woff    :: Vector{Float64}  # 3
    Won     :: Vector{Float64}  # 4
    Solar   :: Vector{Float64}  # 5
    Bio     :: Vector{Float64}  # 6
    BatCap  :: Vector{Float64}  # 7
    BatPow  :: Vector{Float64}  # 8start_year
end

function load_ise_installed_power(hdf5_path)
    # load data matrix from hdf5 file
    uts, M = load_installed_power_from_hdf5(hdf5_path)

#    id = Dict()
#    for (i,pr) in enumerate(prodtypes)
#        id[pr] = i
#        @infoe i, pr
#    end
#    # Unix time stamps are stored as Float64 in D => to Int

#    # convert Unix time stamps to Date objects
#    dates = unix2datetime.(uts)
#    i1 = id["Wind offshore"]
#    i2 = id["Wind onshore"]
#    i3 = id["Solar"]
#    i4 = id["Biomass"]
#    i5 = id["Battery Storage (Capacity)"]
#    i6 = id["Battery Storage (Power)"]
#     
#    DD = Matrix{Float64}(undef, size(D,1), 6)
#    for (i, j) in enumerate((i1, i2, i3, i4, i5, i6))
#        DD[:,i] = D[:,j]
#    end

    uts, M
end

function load_and_iterpolate_installed_power(data_dir, punit, uts_pubpower)
    hp = @sprintf("installed_power.hdf5")
    hdf5_path = joinpath(data_dir, hp)
    uts_instpower, installed_power = load_ise_installed_power(hdf5_path)

    n_uts = length(uts_pubpower)
    n_types = size(installed_power,2)

    #Woff, Won, Solar, Bio, BatCap, BatPow = installed_power

    GW_to_unit = uconversion_factor(u_GW, 1.0*punit)

    installed_power2 = Matrix{Float64}(undef, n_uts, n_types)
    for i in 1:n_types
        interp_linear_extrap = linear_interpolation(uts_instpower, installed_power[:,i], extrapolation_bc=Line())
        installed_power2[:,i] = interp_linear_extrap(uts_pubpower) .* GW_to_unit
    end
    @infoe "n_uts =", n_uts
    @infoe "size(installed_power2) =", size(installed_power2)

    uts_pubpower, installed_power2
end

function InstalledPower(pd::PowerData, data_dir, punit, scale_Bio)
    uts_pubpower, installed_power = load_and_iterpolate_installed_power(data_dir, punit, pd.uts)
    installed_power[:,4] = installed_power[:,4] * scale_Bio

    Woff    = installed_power[:,1]
    Won     = installed_power[:,2]
    Solar   = installed_power[:,3]
    Bio     = installed_power[:,4]
    BatCap  = installed_power[:,5]
    BatPow  = installed_power[:,6]

    A = [Woff, Won, Solar, Bio, BatCap, BatPow]
    B = []
    n_rescaled = 0
    for a in A
        ipmin = minimum(a)
        ipmax = maximum(a)
        if ipmin < 1.0e-10*ipmax
            ipmin = 1.0e-10*ipmax
            push!(B,  @. ifelse(a < ipmin, ipmin, installed_power[:,i]))
            n_rescaled += 1
        else
            push!(B, a)
        end
    end
    @infoe @sprintf("Installed power: # ipmin < 1.0e-10*ipmax = %d", n_rescaled)
    #                      Woff   Won   Sol   Bio   BatCap   BatPow
    InstalledPower(pd.dates, pd.uts, B[1], B[2], B[3], B[4], B[5], B[6])
end
