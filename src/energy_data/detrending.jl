struct DetrendedPowerData
    dates       :: Vector{DateTime} # 1
    uts         :: Vector{Int64}    # 2
    Load        :: Vector{Float64}  # 3
    Woff        :: Vector{Float64}  # 4
    Won         :: Vector{Float64}  # 5 
    Solar       :: Vector{Float64}  # 6 
    Bio         :: Vector{Float64}  # 7 
    WWSBPower   :: Vector{Float64}  # 8 
    Load_trend  :: Vector{Float64}  # 9
    Woff_trend  :: Vector{Float64}  # 10 
    Won_trend   :: Vector{Float64}  # 11 
    Solar_trend :: Vector{Float64}  # 12 
    Bio_trend   :: Vector{Float64}  # 13 
    WWSB_trend  :: Vector{Float64}  # 14 
end

function DetrendedPowerData(df)
    DetrendedPowerData(
        df[!,:dates      ],
        df[!,:uts        ],
        df[!,:Load       ],
        df[!,:Woff       ],
        df[!,:Won        ],
        df[!,:Solar      ],
        df[!,:Bio        ],
        df[!,:WWSBPower  ],
        df[!,:Load_trend ],
        df[!,:Woff_trend ],
        df[!,:Won_trend  ],
        df[!,:Solar_trend],
        df[!,:Bio_trend  ],
        df[!,:WWSB_trend ])
end

"""
    Woff, Won, Solar are scaled such that Woff+Won+Solar+Bio = Load
    Bio assumed not to increase further
"""
function scale_wind_and_solar_by_load!(Load, Woff, Won, Solar, Bio)
    WWS  = Woff .+ Won .+ Solar

    mean_Load = mean(Load)
    mean_Bio  = mean(Bio)
    mean_WWS  = mean(WWS)
    # mean_L*op = mean_R * scale + mean_B
    scale = (mean_Load - mean_Bio) / mean_WWS

    Woff  = Woff  .* scale
    Won   = Won   .* scale
    Solar = Solar .* scale

    Woff, Won, Solar, scale
end

"""
"""
function detrend_time_series(a, trend)
    nh = div(length(a),2)
    a_de = @. a * trend[nh] / trend
    sum_a = sum(a)
    sum_a_de = sum(a_de)
    a_de = @. a_de * sum_a/sum_a_de
    a_de
end

function scale_by_installed_power(public_power, installed_power)
    IP_Woff  = installed_power.Woff  ./ mean(installed_power.Woff)
    IP_Won   = installed_power.Won   ./ mean(installed_power.Won)
    IP_Solar = installed_power.Solar ./ mean(installed_power.Solar)
    IP_Bio   = installed_power.Bio   ./ mean(installed_power.Bio)

    Woff_ip  = @. public_power.Woff  / IP_Woff
    Won_ip   = @. public_power.Won   / IP_Won
    Solar_ip = @. public_power.Solar / IP_Solar
    Bio_ip   = @. public_power.Bio   / IP_Bio
 
    Woff_ip  .*= mean(public_power.Woff)  ./ mean(Woff_ip)
    Won_ip   .*= mean(public_power.Won)   ./ mean(Won_ip)
    Solar_ip .*= mean(public_power.Solar) ./ mean(Solar_ip)
    Bio_ip   .*= mean(public_power.Bio)   ./ mean(Bio_ip)

    Woff_ip, Won_ip, Solar_ip, Bio_ip
end

function get_detrended_public_power(public_power, installed_power, par)
    Woff, Won, Solar, Bio = if par.scale_with_installed_power_p
        scale_by_installed_power(public_power, installed_power)
    else
        public_power.Woff, public_power.Won, public_power.Solar, public_power.Bio  
    end

    uts = public_power.uts

    Load_trend  = polynomial_fit(uts, public_power.Load; degree=2)(uts)
    Load_de = detrend_time_series(public_power.Load, Load_trend)

    Woff_trend  = polynomial_fit(uts, Woff  .* s1; degree=2)(uts)
    Won_trend   = polynomial_fit(uts, Won   .* s2; degree=2)(uts)
    Solar_trend = polynomial_fit(uts, Solar .* s3; degree=2)(uts)
    Bio_trend   = polynomial_fit(uts, Bio   .* s4; degree=2)(uts)

    Woff_de  = detrend_time_series(Woff , Woff_trend )
    Won_de   = detrend_time_series(Won  , Won_trend  )
    Solar_de = detrend_time_series(Solar, Solar_trend)
    Bio_de   = detrend_time_series(Bio  , Bio_trend  )

    Woff_de, Woff_trend, Won_de, Won_trend, Solar_de, Solar_trend, Bio_de, Bio_trend 

    Woff_sc, Won_sc, Solar_sc, scale = scale_wind_and_solar_by_load!(Load_de, Woff_de, Won_de, Solar_de, Bio_de)
    @info @sprintf("scale_factor = %f", scale)

    WWSB = Won_sc .+ Woff_sc .+ Solar_sc .+ Bio_de
    WWSB_trend  = polynomial_fit(uts, WWSB; degree=2)(uts)
    WWSB_de = detrend_time_series(WWSB, WWSB_trend)

    Load_sum = sum(Load_de)
    WWSB_sum = sum(WWSB_de)
    @info @sprintf("Load_sum = %10.4e, WWSB_sum = %10.4e, sum_L-sum_P = %10.4e", Load_sum, WWSB_sum, Load_sum-WWSB_sum)

    DataFrame([
        :dates,
        :uts,
        :Load,
        :Woff,
        :Won,
        :Solar,
        :Bio,
        :WWSBPower,
        :Load_trend  ,
        :Woff_trend  ,
        :Won_trend   ,
        :Solar_trend ,
        :Bio_trend   ,
        :WWSB_trend  ] 
        .=> 
        [public_power.dates,
        uts,
        Load_de,
        Woff_sc,
        Won_sc,
        Solar_sc,
        Bio_de,
        WWSB_de,
        Load_trend  ,
        Woff_trend  ,
        Won_trend   ,
        Solar_trend ,
        Bio_trend   ,
        WWSB_trend  ]
    )
end
