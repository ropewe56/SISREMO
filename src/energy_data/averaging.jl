"""
    hours_to_averaging_steps(dates, averaging_hours)

    returns the number os steps and dates at the staps
"""
function hours_to_averaging_steps(dates, averaging_hours)
    # time difference in ms
    ms = Dates.value(dates[2] - dates[1])
    # time difference in hours
    Δh = ms/(3.6e6)
    c = floor(Int64, 1.0/Δh)
    steps = averaging_hours * c
    steps, dates[1:steps:end]
end

"""
    mean_average(data, dates, averaging_hours)

    average data and dates over averaging_hours
"""
function mean_average(data, dates, averaging_hours)
    averaging_steps, avdates = hours_to_averaging_steps(dates, averaging_hours)

    Nh = div(averaging_steps, 2)
    NE = size(data, 1)

    avdata = []
    avdates = []
    lio = 1
    while lio < NE - averaging_steps
        push!(avdata, mean(data[lio:lio+averaging_steps]))
        push!(avdates, dates[lio + Nh])
        lio = lio + averaging_steps
    end

    avdata, avdates
end

"""
    moving average rectangle window
    y  : input Vector
    iw : number of averaging points
"""
function moving_average_(y::Vector{Float64}, iw::Int64)
    n = length(y)
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
    avdata = moving_average_(data, averaging_steps)
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

function plot_averaged(averaged_power, averaging_hours)
    dates = averaged_power.dates
    Load  = averaged_power.Load
    Won   = averaged_power.Won
    Woff  = averaged_power.Woff
    Solar = averaged_power.Solar
    WWSB  = Won .+ Woff .+ Solar .+ averaged_power.Bio

    plt.figure()
    plt.plot(dates, Load, label="Load")
    plt.plot(dates, WWSB, label="WWSB")
    plt.plot(dates, Won .+ Woff, label="Wind")
    plt.plot(dates, Solar, label="Solar")
    plt.plot(dates, WWSB, label="WWSB")
    plt.xlabel("time")
    plt.ylabel("P [GW]")
    plt.legend()
    plt.title(@sprintf("EE averages over %d hours", averaging_hours))
    #plt.savefig(joinpath(FIGDIR, "EE_averaged.png"))

    plt.figure()
    plt.plot(dates, WWSB ./ Load, label="WWSB/Load")
    plt.xlabel("time")
    plt.ylabel("P [GW]")
    plt.legend()
end
