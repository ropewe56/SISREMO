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
    dates, uts, D[:,i1], D[:,i2], D[:,i3], D[:,i4], D[:,i5],  D[:,i6]
end

struct PowerData
    dates     :: Vector{DateTime} # 1
    uts       :: Vector{Int64}    # 2
    Load      :: Vector{Float64}  # 3
    Woff      :: Vector{Float64}  # 4
    Won       :: Vector{Float64}  # 5
    Solar     :: Vector{Float64}  # 6
    Bio       :: Vector{Float64}  # 7
    Nuclear   :: Vector{Float64}  # 8
    WWSBPower :: Vector{Float64}  # 8
end

"""
    EnergyData(eunit)

    unix_seconds
    Biomass
    Cross_border_electricity_trading
    Fossil_brown_coal_lignite
    Fossil_coal_derived_gas
    Fossil_gas
    Fossil_hard_coal
    Fossil_oil
    Geothermal
    Hydro_Run_of_River
    Hydro_pumped_storage
    Hydro_pumped_storage_consumption
    Hydro_water_reservoir
    Load
    Nuclear
    Others
    Renewable_share_of_generation
    Renewable_share_of_load
    Residual_load
    Solar
    Waste
    Wind_offshore
    Wind_onshore

    EnergyData constructor
"""
function PowerData(start_year, end_year, par)
    total_power = get_total_power(start_year, end_year)
    for n in names(total_power) println(n)  end

    Geo_4     = total_power[!,:Geothermal]
    Load_4    = total_power[!,:Load]
    Woff_4    = total_power[!,:Wind_offshore]
    Won_4     = total_power[!,:Wind_onshore]
    Solar_4   = total_power[!,:Solar]
    Bio_4     = total_power[!,:Biomass]
    Nuclear_4 = total_power[!,:Nuclear]

    MW_to_unit = uconversion_factor(par.punit, 1u_MW)

    Geo     = average_to_hour(Geo_4)      .* MW_to_unit .* par.Geo_scale
    Load    = average_to_hour(Load_4)     .* MW_to_unit .* par.Load_scale
    Woff    = average_to_hour(Woff_4)     .* MW_to_unit .* par.Woff_scale
    Won     = average_to_hour(Won_4)      .* MW_to_unit .* par.Won_scale
    Solar   = average_to_hour(Solar_4)    .* MW_to_unit .* par.Solar_scale
    Bio     = average_to_hour(Bio_4)      .* MW_to_unit .* par.Bio_scale
    Nuclear = average_to_hour(Nuclear_4)  .* MW_to_unit

    WWSBPower = Woff .+ Won .+ Solar .+ Bio

    nl  = length(Load)
    uts = uts_4[1:4:end][1:nl]
    dates = Dates.unix2datetime.(uts)

    @infoe @sprintf("# timesteps = %d, length(Load) = %d, energy conversion = %e", length(uts), length(Load), MW_to_unit)

    PowerData(dates, uts, Load, Woff, Won, Solar, Bio, Nuclear, WWSBPower)
end

struct InstalledPowerData
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
    uts, prodtypes = load_installed_power_from_hdf5(hdf5_path)
    uts, prodtypes
end

function load_and_iterpolate_installed_power(hdf5_dir, punit, uts_pubpower, start_year, end_year)
    total_power = get_total_power(start_year, end_year)
    installed_power = get_installed_power(start_year, end_year)
    uts = [datetime2unix(DateTime(x, 7, 1)) for x in installed_power[!, :time]]

    uts_pubpower = total_power[!,:unix_seconds]

    cols = """time, Battery_storage_capacity, Battery_storage_power, Biomass, Fossil_brown_coal_lignite, 
        Fossil_gas, Fossil_hard_coal, Fossil_oil, Hydro, Hydro_pumped_storage, Nuclear, Other_non_renewable, 
        Solar_AC, Solar_DC, Solar_planned_EEG_2023, Wind_offshore, Wind_onshore, Wind_onshore_planned_EEG_2023"""
    syms = [Symbol(strip(a)) for a in split(cols, ",")]

    GW_to_unit = uconversion_factor(u_GW, 1.0*punit)
    installed_power2 = Dict()

    for sym in syms
        interp_linear_extrap = linear_interpolation(uts, installed_power[!,sym], extrapolation_bc=Line())
        installed_power2[sym]  = interp_linear_extrap(uts_pubpower) .* GW_to_unit
    end
    @infoe "n_uts =", length(uts_pubpower)
    @infoe "size(installed_power2) =", length(keys(installed_power2))

    uts_pubpower, installed_power2
end

function InstalledPowerData(hdf5_dir, power_data, par)
    uts_pubpower, installed_power = load_and_iterpolate_installed_power(hdf5_dir, par.punit, power_data.uts)

    name_map = Dict("Bio"    => "Biomass",
                    "Woff"   => "Wind onshore",
                    "Won"    => "Wind offshore",
                    "Solar"  => "Solar",
                    "BatPow" => "Battery Storage (Power)",
                    "BatCap" => "Battery Storage (Capacity)")
    
    Woff    = installed_power[name_map["Woff"  ]] * par.Woff_scale
    Won     = installed_power[name_map["Won"   ]] * par.Won_scale
    Solar   = installed_power[name_map["Solar" ]] * par.Solar_scale
    Bio     = installed_power[name_map["Bio"   ]] * par.Bio_scale
    BatCap  = installed_power[name_map["BatCap"]]
    BatPow  = installed_power[name_map["BatPow"]]

    A = [Woff, Won, Solar, Bio, BatCap, BatPow]
    B = []
    n_rescaled = 0
    for a in A
        ipmin = minimum(a)
        ipmax = maximum(a)
        if ipmin < 1.0e-10*ipmax
            ipmin = 1.0e-10*ipmax
            a2 = @. ifelse(a < ipmin, ipmin, a)
            push!(B,  a2)
            n_rescaled += 1
        else
            push!(B, a)
        end
    end
    @infoe @sprintf("Installed power: # ipmin < 1.0e-10*ipmax = %d", n_rescaled)
    @infoe @sprintf("interpolated installed data %d", length(B[1]))
    #                                                    Woff  Won   Sol  Bio   BatCap BatPow
    InstalledPowerData(power_data.dates, power_data.uts, B[1], B[2], B[3], B[4], B[5], B[6])
end


