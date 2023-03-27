using CurveFit
using Dates
using JSON

#macro infoe(message)
#    fl = string(basename(string(__source__.file)), ":", (__source__.line))
#    return :( @printf("[ Info: %s | %s\n", $fl, $(esc(message))) )
#end

#linspace(a, b, n) = collect(LinRange(a, b, n))

#"""
#    unix utstamp vergangene Sekunden seit Donnerstag, 1. Januar 1970, 00:00 Uhr UTC.
#    date -u -d @1234567890
#    "2021-01-01T17:00Z" >= T-time=17:00, Z=UTC
#"""
#uts2date(uts)  = @. unix2datetime(uts)
#date2uts(date) = @. floor(Int64, datetime2unix(date))
#
#date2ymwdhms(dt) = year(dt), month(dt), week(dt), day(dt), hour(dt), minute(dt), second(dt)
#
#"""
#    JSON write
#"""
#function to_json(path::String, d)
#    open(path, "w") do out
#        JSON.print(out, d, 4)
#    end
#end
#"""
#    JSON read
#"""
#function from_json2(path)
#    open(path, "r") do ein
#        d = read(ein, String)
#        return JSON.parse(d)
#    end
#end


"""
    polynomial_fit(y, k)
    y : data to fit
    k : order of polynomial
"""
function polynomial_fit(y, k)
    N = length(y)
    x = linspace(0.0, Float64(N), N)
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

function averaging(data, dates, averaging_hours; method = "moving_average")
    if method == "moving_average"
        moving_average(data, dates, averaging_hours)
    elseif method == "mean"
        mean_averaging(data, dates, averaging_hours)
    else
        @infoe @sprintf("method %s not vaialable", method)
    end
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

    Load_avmo = moving_average(Load, averaging_hours)
    EE_avmo   = moving_average(EE, averaging_hours)

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

    @infoe @sprintf("sum(diff2) = %e", sum(diff2))
    @infoe @sprintf("sum(diff3) = %e", sum(diff3))

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
