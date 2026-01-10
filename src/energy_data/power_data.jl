function uconversion_factor(eto, efrom)
    @info eto, efrom
    Float64(getproperty(uconvert(eto, efrom), :val))
end

function save_to_arrow(df::DataFrame, path_arrow)
    Arrow.write(path_arrow, df)
end

function load_from_arrow(path_arrow)
    copy(Arrow.Table(path_arrow) |> DataFrame)
end

function average_to_hour(X)
    n = length(X)
    @. (X[1:4:n-3] + X[2:4:n-2] + X[3:4:n-1] + X[4:4:n]) * 0.25
end

Base.@kwdef mutable struct PowerParameter
    punit                        = u_GW
    scale_Bio                    = 1.0
    Load_scale                   = 1.0
    Woff_scale                   = 1.0
    Won_scale                    = 1.0
    Solar_scale                  = 1.0
    Bio_scale                    = 1.0
    Geo_scale                    = 1.0
    SF1_factor                   = 1.0
    scale_to_installed_power_p   = true
    scale_with_installed_power_p = false
    averaging_hours              = 24*7
    averaging_method             = [:moving_average, :mean][1]
    plot_p                       = true
    plot_all_p                   = false
    log_p                        = false
    sisremo_root                 = SISREMOROOT
    data_dir                     = DATAROOT
    fig_dir                      = FIGDIR
    start_year                   = 2016
    end_year                     = 2025
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

    EnergyData constructor
    par = PowerParameter()
    plt.plot(pp[!,:unix_seconds], pp[!,:Solar])
    plt.plot(uts, average_to_hour(pp[!,:Solar]))
"""
function get_public_power_data(date1, date2, par)
    pp = select_from_db(date1, date2, ["public_power"])["public_power"]

    # energy charts data are in MW, MW_to_unit is conversion factor to eunit (MW, GW, TW)
    MW_to_unit = uconversion_factor(par.punit, 1u_MW)

    df = DataFrame()
    df[!,:uts] = average_to_hour(pp[!,:unix_seconds])
    df[!,:dates] = Dates.unix2datetime.(df[!,:uts])

    df[!,:Load   ] = average_to_hour(pp[!,:Load])          .* MW_to_unit .* par.Load_scale
    df[!,:Woff   ] = average_to_hour(pp[!,:Wind_offshore]) .* MW_to_unit .* par.Woff_scale
    df[!,:Won    ] = average_to_hour(pp[!,:Wind_onshore])  .* MW_to_unit .* par.Won_scale
    df[!,:Solar  ] = average_to_hour(pp[!,:Solar])         .* MW_to_unit .* par.Solar_scale
    df[!,:Bio    ] = average_to_hour(pp[!,:Biomass])       .* MW_to_unit .* par.Bio_scale
    df[!,:Nuclear] = average_to_hour(pp[!,:Nuclear])       .* MW_to_unit

    df[!,:WWSBPower] = df[!,:Woff] .+ df[!,:Won] .+ df[!,:Solar] .+ df[!,:Bio]

    df
end

function PowerData(df)
    PowerData(
        df[!,:dates    ],
        df[!,:uts      ],
        df[!,:Load     ],
        df[!,:Woff     ],
        df[!,:Won      ],
        df[!,:Solar    ],
        df[!,:Bio      ],
        df[!,:Nuclear  ],
        df[!,:WWSBPower])
end

struct InstalledPowerData
    dates   :: Vector{DateTime} # 1
    uts     :: Vector{Int64}    # 2
    Woff    :: Vector{Float64}  # 3
    Won     :: Vector{Float64}  # 4
    Solar   :: Vector{Float64}  # 5
    Bio     :: Vector{Float64}  # 6
    Nuclear :: Vector{Float64}  # 6
    WWSBPower:: Vector{Float64}  # 6
    BatCap  :: Vector{Float64}  # 7
    BatPow  :: Vector{Float64}  # 8start_year
end

"""
    uts_pubpower = power_data.uts
    punit = par.punit
"""
function iterpolate_installed_power_1(installed_power_y::DataFrame, punit, uts_pubpower)
    year = installed_power_y[!,:time]
    uts = [datetime2unix(DateTime(@sprintf("%d-07-01", t))) for t in year]

    values = Vector{Any}([uts_pubpower])
    colnames = names(installed_power_y)
    for k in colnames[2:end]
        intp = linear_interpolation(uts, installed_power_y[!,k], extrapolation_bc=Line())
        push!(values, intp(uts_pubpower))
    end
    df = DataFrame(colnames .=> values)
    df
end

"""
    uts_pubpower = power_data.uts
    punit = par.punit
    uts_pubpower = power_data.uts
"""
function iterpolate_installed_power_2(installed_power_y::DataFrame, uts_pubpower)
    yip = installed_power_y[!,:time]
    colnames = names(installed_power_y)
    nc = ncol(installed_power_y)

    yp = Dates.year.(unix2datetime.(uts_pubpower))

    dfi = DataFrame(["i", "y"] .=> [1:length(yp), yp])

    values = Vector{Any}(undef, 0)
    push!(values, uts_pubpower)
    for i in 2:nc
        push!(values, zeros(Float64, length(uts_pubpower)))
    end
    dfIP = DataFrame(names(installed_power_y) .=> values)

    for (iy, y) in enumerate(yip)
        ii = filter(row -> row[:y] == y, dfi)[!,:i]
        if length(ii) > 0
            u = collect(range(0.0, 1.0, length(ii)))
            for k in colnames[2:end]       
                dfIP[ii[1]:ii[end], k] = @. (1.0-u) * installed_power_y[iy,k] + u * installed_power_y[iy+1,k]
            end
        end
    end
    dfIP
end

function get_installed_power_data(power_data, par)
    date1, date2 = DateTime("2016-01-01"), now()
    installed_power_y = select_from_db(date1, date2, ["installed_power"])["installed_power"]

    installed_power_h = iterpolate_installed_power_2(installed_power_y, power_data.uts)

    GW_to_unit = uconversion_factor(u_GW, 1.0*par.punit)
    df = DataFrame()

    df[!,:uts]    = installed_power_h[!,:time]
    df[!,:dates]  = unix2datetime.(df[!,:uts])

    df[!,:Woff]      = installed_power_h[!,:Wind_offshore] * par.Woff_scale  * GW_to_unit
    df[!,:Won]       = installed_power_h[!,:Wind_onshore]  * par.Won_scale   * GW_to_unit
    df[!,:Solar]     = (installed_power_h[!,:Solar_AC] .+ installed_power_h[!,:Solar_DC]) .* (par.Solar_scale * GW_to_unit)
    df[!,:Bio]       = installed_power_h[!,:Biomass] * par.Bio_scale   * GW_to_unit
    df[!,:Nuclear]   = installed_power_h[!,:Nuclear]  * GW_to_unit
    df[!,:WWSBPower] = df[!,:Woff] .+ df[!,:Won] .+  df[!,:Solar] .+ df[!,:Bio]
    df[!,:BatCap]    = installed_power_h[!,:Battery_storage_capacity] * GW_to_unit
    df[!,:BatPow]    = installed_power_h[!,:Battery_storage_power] * GW_to_unit

    df
end

function InstalledPowerData(df)
    InstalledPowerData(
        df[!,:dates ],
        df[!,:uts   ],
        df[!,:Woff  ],
        df[!,:Won   ],
        df[!,:Solar ],
        df[!,:Bio   ],
        df[!,:Nuclear],
        df[!,:WWSBPower],
        df[!,:BatCap],
        df[!,:BatPow])
end

struct AveragedPowerData
    dates    :: Vector{DateTime} # 1
    uts      :: Vector{Int64}    # 2
    Load     :: Vector{Float64}  # 3
    Woff     :: Vector{Float64}  # 4
    Won      :: Vector{Float64}  # 5
    Solar    :: Vector{Float64}  # 6
    Bio      :: Vector{Float64}  # 7
    Nuclear  :: Vector{Float64}  # 8
    WWSBPower:: Vector{Float64}  # 9 # sum of Woff Won Solar Bio
end

function get_averaged_power_data(power_data, averaging_hours, averaging_method)
    Load     , dates_av = averaging(power_data.Load     , power_data.dates, averaging_hours, method = averaging_method)
    Woff     , dates_av = averaging(power_data.Woff     , power_data.dates, averaging_hours, method = averaging_method)
    Won      , dates_av = averaging(power_data.Won      , power_data.dates, averaging_hours, method = averaging_method)
    Solar    , dates_av = averaging(power_data.Solar    , power_data.dates, averaging_hours, method = averaging_method)
    Bio      , dates_av = averaging(power_data.Bio      , power_data.dates, averaging_hours, method = averaging_method)
    Nuclear  , dates_av = averaging(power_data.Nuclear  , power_data.dates, averaging_hours, method = averaging_method)
    WWSBPower, dates_av = averaging(power_data.WWSBPower, power_data.dates, averaging_hours, method = averaging_method)

    
    uts_av = [Dates.datetime2unix(x) for x in dates_av]

    DataFrame([:dates, :uts, :Load, :Woff, :Won, :Solar, :Bio, :Nuclear, :WWSBpower] .=>
                [dates_av, uts_av, Load, Woff, Won, Solar, Bio, Nuclear, WWSBPower])
end

function AveragedPowerData(df)
    AveragedPowerData(
        df[!,:dates    ],
        df[!,:uts      ],
        df[!,:Load     ],
        df[!,:Woff     ],
        df[!,:Won      ],
        df[!,:Solar    ],
        df[!,:Bio      ],
        df[!,:Nuclear  ],
        df[!,:WWSBPower])
end