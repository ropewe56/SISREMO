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
function get_public_public_power(date1, date2, par)
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
    uts = public_power.uts
    punit = par.punit
"""
function iterpolate_installed_power_2(installed_power::DataFrame, uts)
    yip = installed_power[!,:time]
    utsip = datetime2unix.(DateTime.(yip, 1, 1))

    dates = unix2datetime.(uts)

    dfip = DataFrame()
    dfip[!,:dates] = dates
    dfip[!,:uts]  = uts
    dfip[!,:Won]  = polynomial_fit(utsip, installed_power[!,:Wind_onshore]; degree = 2)(uts)
    dfip[!,:Woff] = polynomial_fit(utsip, installed_power[!,:Wind_offshore]; degree = 2)(uts)
    SAC = polynomial_fit(utsip, installed_power[!,:Solar_AC]; degree = 2)(uts)
    SDC = polynomial_fit(utsip, installed_power[!,:Solar_DC]; degree = 2)(uts)
    dfip[!,:Solar] = SAC .+ SDC
    dfip[!,:Bio]  = polynomial_fit(utsip, installed_power[!,:Biomass]; degree = 2)(uts)

    dfip[!,:Nuclear] = polynomial_fit(utsip, installed_power[!,:Nuclear]; degree = 2)(uts)
    dfip[!,:WWSBPower] = dfip[!,:Won] .+ dfip[!,:Woff] .+   dfip[!,:Bio] .+ dfip[!,:Solar] 

    dfip[!,:BatCap] = polynomial_fit(utsip, installed_power[!,:Battery_storage_capacity]; degree = 2)(uts)
    dfip[!,:BatPow] = polynomial_fit(utsip, installed_power[!,:Battery_storage_power]; degree = 2)(uts)

    dfip
end

function get_installed_public_power(public_power, par)
    date1, date2 = DateTime("2016-01-01"), now()
    installed_power_y = select_from_db(date1, date2, ["installed_power"])["installed_power"]

    dfip = iterpolate_installed_power_2(installed_power_y, public_power.uts)

    GW_to_unit = uconversion_factor(u_GW, 1.0*par.punit)
    dfip[!,:Woff]      = dfip[!,:Woff]      .* (par.Woff_scale  * GW_to_unit)
    dfip[!,:Won]       = dfip[!,:Won]       .* (par.Won_scale   * GW_to_unit)
    dfip[!,:Solar]     = dfip[!,:Solar]     .* (par.Solar_scale * GW_to_unit)
    dfip[!,:Bio]       = dfip[!,:Bio]       .* (par.Bio_scale   * GW_to_unit)
    dfip[!,:Nuclear]   = dfip[!,:Nuclear]   .* GW_to_unit
    dfip[!,:WWSBPower] = dfip[!,:WWSBPower] .* GW_to_unit
    dfip[!,:BatCap]    = dfip[!,:BatCap]    .* GW_to_unit
    dfip[!,:BatPow]    = dfip[!,:BatPow]    .* GW_to_unit

    dfip
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

function get_averaged_public_power(public_power, averaging_hours, averaging_method)
    Load     , dates_av = averaging(public_power.Load     , public_power.dates, averaging_hours, method = averaging_method)
    Woff     , dates_av = averaging(public_power.Woff     , public_power.dates, averaging_hours, method = averaging_method)
    Won      , dates_av = averaging(public_power.Won      , public_power.dates, averaging_hours, method = averaging_method)
    Solar    , dates_av = averaging(public_power.Solar    , public_power.dates, averaging_hours, method = averaging_method)
    Bio      , dates_av = averaging(public_power.Bio      , public_power.dates, averaging_hours, method = averaging_method)
    Nuclear  , dates_av = averaging(public_power.Nuclear  , public_power.dates, averaging_hours, method = averaging_method)
    WWSBPower, dates_av = averaging(public_power.WWSBPower, public_power.dates, averaging_hours, method = averaging_method)

    
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