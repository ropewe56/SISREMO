using CurveFit
using Dates

"""
    polynomial_fit(y, k)
    y : data to fit
    k : order of polynomial
"""
function polynomial_fit(y, k)
    N = length(y)
    x = collect(range(0.0, Float64(N), N))
    p = poly_fit(x, y, k)
    yf = p[1]
    for i in 1:k
        yf = yf .+ p[i+1] .* x.^i
    end
    yf
end

function spline_approximation(y)
    N = length(y)
    x = linspace(0.0, Float64(N), N)
    #cs = scp.CubicSpline(x,y)
    cs = scp.UnivariateSpline(x, y, k=3, s=10.0, ext=0, check_finite=false)
    y = cs(x)
    pl.figure()
    pl.plot(y, label="scipy")
    pl.legend()
    y
end

function hours_to_averaging_steps(dates, averaging_hours)
    ms = Dates.value(dates[2] - dates[1])
    Δh = ms/(3.6e6)
    c = floor(Int64, 1.0/Δh)
    steps = averaging_hours * c
    steps, dates[1:steps:end]
end

function mean_averaging(data, dates, averaging_hours)
    averaging_steps, avdates = hours_to_averaging_steps(dates, averaging_hours)
    Nh = div(averaging_steps, 2)
    NE = size(data, 1)
    lio = 1

    avdata = []
    avdates = []
    while lio < NE - averaging_steps
        push!(avdata, mean(data[lio:lio+averaging_steps]))
        push!(avdates, dates[lio + Nh])
        lio = lio + averaging_steps
    end

    avdata, avdates
end

"""
    moving average rect window
    y  : input Vector
    iw : number of averaging points
"""
function moving_averaging(y::Vector{Float64}, iw::Int64)
    n = size(y,1)
    ma = zeros(Float64, n)
    mi = zeros(Float64, n)
    ma[1] = sum(y[1:iw])
    mi[1] = Float64(iw)
    ym  = 0.0
    mam = 0.0
    mim = 0.0

    for i in 2:n
        im = i - iw
        ip = i + iw

        dma1 = 0.0
        dmi1 = 0.0
        dma2 = 0.0
        dmi2 = 0.0

        if im >= 1
            dma1 = - y[im]
            dmi1 = - 1.0
        end
        if ip <= n
            dma2 = dma2 + y[ip]
            dmi2 = dmi2 + 1.0
        end
        ma[i] = ma[i-1] + dma1 + dma2
        mi[i] = mi[i-1] + dmi1 + dmi2

        ym += y[i]
    end

    for i in 1:n
        if mi[i] > 0.0
            ma[i] = ma[i] / mi[i]
        end
        mam += ma[i]
        mim += mi[i]
    end
    ma
end

function moving_average(data, dates, averaging_hours)
    averaging_steps, avdates = hours_to_averaging_steps(dates, averaging_hours)
    avdata = moving_averaging(data, averaging_steps)
    avdata, avdates
end

function averaging(data, dates, averaging_hours; method = :moving_average)
    avdata, avdates = if method == :moving_average
        moving_average(data, dates, averaging_hours)
    elseif method == :mean
        mean_averaging(data, dates, averaging_hours)
    else
        @info @sprintf("method %s not vaialable", method)
    end

    if length(avdata) == length(avdates)
        return avdata, avdates
    elseif length(avdata) == length(dates)
        return avdata, dates
    end
    @warn @sprintf("length of data and dates don't fit.")
    return avdata, avdates
end

function simple_damping(y, α)
    n = size(y,1)
    L = Vector{Float64}(undef, n)
    L[1] = y[1]
    for i in 2:n
        L[i] = α*y[i] + (1.0-α)*L[i-1]
    end
    L
end

function average_data(Load, EE, dates, averaging_hours)
    Load_avme, dates_avme = averaging_mean(Load, dates, averaging_hours)
    EE_avme,   dates_avme = averaging_mean(EE, dates, averaging_hours)

    Load_avmo = moving_average(Load, dates_avme, averaging_hours)
    EE_avmo   = moving_average(EE, dates_avme, averaging_hours)

    Load_avme, Load_avmo, EE_avme, EE_avmo, dates_avme
end

function plot_averaged(Load_avme, Load_avmo, EE_avme, EE_avmo, dates_avme, dates, averaging_hours, fig_dir)
    pl.figure()
    pl.plot(dates,      EE_avmo, label="EE_avmo")
    pl.plot(dates_avme, EE_avme, label="EE_avme")
    pl.xlabel("time")
    pl.ylabel("P [MW]")
    pl.legend()
    pl.title(@sprintf("EE averages over %d hours", averaging_hours))
    pl.savefig(joinpath(fig_dir, "EE_averaged.png"))

    pl.figure()
    pl.plot(dates,      Load_avmo, label="Load_avmo")
    pl.plot(dates_avme, Load_avme, label="Load_avme")
    pl.xlabel("time")
    pl.ylabel("P [MW]")
    pl.legend()
    pl.title(@sprintf("Load averages over %d hours", averaging_hours))
    pl.savefig(joinpath(fig_dir, "Load_averaged.png"))
end

# Detrend time series data
# * https://en.wikipedia.org/wiki/Singular_spectrum_analysis
# * https://github.com/baggepinnen/SingularSpectrumAnalysis.jl
# window_length = 7*24
# yt, ys = SingularSpectrumAnalysis.analyze(EEmc, window_length) # trend and seasonal components
function linfit(A)
    N   = length(A)
    t   = linspace(0.0, Float64(N), N)
    tm  = mean(t)
    t2m = mean(t.^2)
    Am  = mean(A)
    Atm = mean(A.*t)

    B = Matrix{Float64}(undef, 2, 2)
    B[1,1] = 1.0
    B[2,1] = tm
    B[1,2] = tm
    B[2,2] = t2m

    y = Vector{Float64}(undef, 2)
    y[1] = Am
    y[2] = Atm

    x = linsys22(B, y)
    a = x[1]
    b = x[2]

    a, b, t
end

function detrend_EE(Load, EE, Load_avme, Load_avmo, EE_avme, EE_avmo, dates_avme, dates, fig_dir)

    EE_avmo_scaled = EE_avmo * mean(Load)/mean(EE)

    EE_avmo_scaled_mean = mean(EE_avmo_scaled)

    EE_trend2 = polynomial_fit(EE_avmo_scaled, 2)
    EE_trend3 = polynomial_fit(EE_avmo_scaled, 3)

    EE_avma_scaled_de2 = EE_avmo_scaled .- EE_trend2 .+ EE_avmo_scaled_mean
    EE_avma_scaled_de3 = EE_avmo_scaled .- EE_trend3 .+ EE_avmo_scaled_mean

    ls = sum(Load_avmo)

    es2 = sum(EE_avma_scaled_de2)
    EE_avma_scaled_de2 = EE_avma_scaled_de2 .* (ls/es2)
    diff2 = (EE_avma_scaled_de2 - Load_avmo) # MW

    es3 = sum(EE_avma_scaled_de3)
    EE_avma_scaled_de3 = EE_avma_scaled_de3 .* (ls/es3)
    diff3 = (EE_avma_scaled_de3 - Load_avmo) # MW

    @info @sprintf("sum(diff2) = %e", sum(diff2))
    @info @sprintf("sum(diff3) = %e", sum(diff3))

    pl.figure()
    pl.plot(dates, 1.0e-6.*EE_avmo_scaled,     label = "EE_avmo_scaled")
    pl.plot(dates, 1.0e-6.*EE_avma_scaled_de2, label = "EE_avma_scaled_de2")
    pl.plot(dates, 1.0e-6.*EE_avma_scaled_de3, label = "EE_avma_scaled_de3")
    pl.plot(dates, 1.0e-6.*EE_trend2,      label = "EE_trend2")
    pl.plot(dates, 1.0e-6.*EE_trend3,      label = "EE_trend3")
    pl.xlabel("time")
    pl.ylabel("P [TW]")
    pl.legend()
    pl.title("detrended")
    pl.savefig(joinpath(fig_dir, "EE_detrended.png"))


    pl.figure()
    pl.plot(dates, 1.0e-6.*diff2, label="EE_avma_scaled_de - Load_avmo 2")
    pl.plot(dates, 1.0e-6.*diff3, label="EE_avma_scaled_de - Load_avmo 3")
    pl.xlabel("time")
    pl.ylabel("(P_E - P_L) [TW]")
    pl.legend()
    pl.grid()
    pl.title("diff detrended")
    pl.savefig(joinpath(fig_dir, "EE_diff_detrended.png"))

    # scaled detrnde 2'nd order fit
    EE_avma_scaled_de2
end

function detrend_time_series(A; pol_order=2)
    A_trend = polynomial_fit(A, pol_order)

    #A_de = @. A * mean(A_trend) / A_trend
    nh = div(length(A),2)
    A_de = @. A * A_trend[nh] / A_trend

    sumA = sum(A)
    sumA_de = sum(A_de)
    A_de = @. A_de * sumA/sumA_de

    A_de, A_trend
end

"""
    Woff, Won, Solar are scaled such that Woff+Won+Solar+Bio = Load
    Bio assumed not to increase further
"""
function scale_power(Load, Woff, Won, Solar, Bio)
    WWS  = Woff .+ Won .+ Solar

    mean_Load = mean(Load)
    mean_Bio  = mean(Bio)
    mean_WWS  = mean(WWS)

    # mean_L*op = mean_R * scale + mean_B
    scale = (mean_Load - mean_Bio) / mean_WWS

    Woff_sc  = Woff  .* scale
    Won_sc   = Won   .* scale
    Solar_sc = Solar .* scale
    WWSB_sc  = Woff_sc .+ Won_sc .+ Solar_sc  .+ Bio

    Woff_sc, Won_sc, Solar_sc, WWSB_sc, scale
end

function renewables_detrend_and_scale(par, power_data)
    local Woff_de 
    local Won_de  
    local Solar_de
    local Bio_de  
    local Load_de 
    local Woff_trend 
    local Won_trend  
    local Solar_trend
    local Bio_trend  
    local Load_trend 
    
    if par.scale_with_installed_power_p
        IP = InstalledPowerData(par, power_data)

        IP_Woff  = @. IP.Woff  / mean(IP.Woff)
        IP_Won   = @. IP.Won   / mean(IP.Won)
        IP_Solar = @. IP.Solar / mean(IP.Solar)
        IP_Bio   = @. IP.Bio   / mean(IP.Bio)

        Woff_ip  = @. power_data.Woff  / IP_Woff
        Won_ip   = @. power_data.Won   / IP_Won
        Solar_ip = @. power_data.Solar / IP_Solar
        Bio_ip   = @. power_data.Bio   / IP_Bio

        s1 = mean(power_data.Woff) / mean(Woff_ip)
        s2 = mean(power_data.Won)  / mean(Won_ip)
        s3 = mean(power_data.Solar)/ mean(Solar_ip)
        s4 = mean(power_data.Bio)  / mean(Bio_ip)

        Woff_de , Woff_trend  = detrend_time_series(Woff_ip  .* s1)
        Won_de  , Won_trend   = detrend_time_series(Won_ip   .* s2)
        Solar_de, Solar_trend = detrend_time_series(Solar_ip .* s3)
        Bio_de  , Bio_trend   = detrend_time_series(Bio_ip   .* s4)
        Load_de , Load_trend  = detrend_time_series(power_data.Load)
    else
        Woff_de , Woff_trend  = detrend_time_series(power_data.Woff)
        Won_de  , Won_trend   = detrend_time_series(power_data.Won)
        Solar_de, Solar_trend = detrend_time_series(power_data.Solar)
        Bio_de  , Bio_trend   = detrend_time_series(power_data.Bio)
        Load_de , Load_trend  = detrend_time_series(power_data.Load)
    end

    Woff_sc, Won_sc, Solar_sc, WWSB_sc, sc = scale_power(Load_de, Woff_de, Won_de, Solar_de, Bio_de)
    @info @sprintf("scale_factor = %f", sc)

    WWSB_de, WWSB_trend = detrend_time_series(WWSB_sc)

    Load_sum = sum(Load_de)
    WWSB_sum = sum(WWSB_sc)
    @info @sprintf("Load_sum = %10.4e, WWSB_sum = %10.4e, sum_L-sum_P = %10.4e", Load_sum, WWSB_sum, Load_sum-WWSB_sum)

    powers_de = DetrendedPowerData(power_data.dates, power_data.uts, 
                    Load_de, Woff_sc, Won_sc, Solar_sc, Bio_de, WWSB_de, 
                    Load_trend, Woff_trend, Won_trend, Solar_trend, Bio_trend, WWSB_trend)
    powers_de
end
