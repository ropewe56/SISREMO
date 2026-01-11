"""
    plot_powers(dates, Load, RP, averaging_hours, fig_dir, punit, fig)

    plot Load and RP

    dates - times
    Load -
    RP   -
    averaging_hours - if data are averaged number of hours to avergae over
    fig_dir - directory where figures are saved to
    punit - (MW, GW, TW)
    fig - number of matplotlib figure
"""
function plot_powers(dates::Vector{DateTime}, Load_ec::Vector{Float64}, Load::Vector{Float64}, 
                        WWSB_ec::Vector{Float64}, WWSB::Vector{Float64},
                        averaging_hours::Int64, fig_dir::String, punit, fig::Vector{Int64})

    mkpath(fig_dir)

    plt.figure(fig[1]); fig[1] += 1
    plt.plot(dates, Load_ec, label="Load")
    plt.plot(dates, WWSB_ec, label="WWSB")
    plt.xlabel("time")
    plt.ylabel(@sprintf("P [%s]", punit))
    plt.grid()
    plt.title("Energy-Charts data")
    plt.legend()

    plt.figure(fig[1]); fig[1] += 1
    plt.plot(dates, WWSB, label="WWSB")
    plt.plot(dates, Load, label="Load")
    plt.xlabel("time")
    plt.ylabel(@sprintf("P [%s]", punit))
    plt.grid()
    plt.title("scaled  data")
    plt.legend()

    if averaging_hours > 1
        plt.title(@sprintf("averaged, %d days window", div(averaging_hours,24)))
        plt.savefig(joinpath(fig_dir, @sprintf("WWSB_averaged.png")))
    else
        @info fig_dir
        plt.savefig(joinpath(fig_dir, @sprintf("WWSB.png")))
    end
end

function plot_detrended(dates::Vector{DateTime}, WWSB::Vector{Float64}, WWSB_de::Vector{Float64}, 
                        ΔEL::Vector{Float64},
                        Load::Vector{Float64}, Load_de::Vector{Float64}, Load_trend::Vector{Float64},
                        punit, fig_dir::String, fig::Vector{Int64}; data_are_averaged = false)
    mkpath(fig_dir)

    label1 = "Load"
    label2 = "Load_de"
    path1  = "Load_detrended.png"
    path2  = "Load_diff_detrended.png"
    if data_are_averaged
        label1 = "Load_av"
        label2 = "Load_av_de"
        path1  = "Load_av_detrended.png"
    end

    plt.figure(fig[1]); fig[1] += 1
    plt.plot(dates, Load      , label = label1)
    plt.plot(dates, Load_de   , label = label2)
    plt.plot(dates, Load_trend, "r", linewidth=3, label = "Load_trend")
    plt.xlabel("time")
    plt.ylabel(@sprintf("P [%s]", punit))
    plt.grid()
    plt.legend()
    plt.title("Load detrended")
    plt.savefig(joinpath(fig_dir, path1))

    label1 = "WWSB"
    label2 = "WWSB_de"
    label3 = "WWSB_de - Load"
    path1  = "WWSB_detrended.png"
    path2  = "WWSB_diff_detrended.png"
    if data_are_averaged
        label1 = "WWSB_av"
        label2 = "WWSB_av_de"
        label3 = "WWSB_av_de - Load_av"
        path1  = "WWSB_av_detrended.png"
        path2  = "WWSB_av_diff_detrended.png"
    end

    plt.figure(fig[1]); fig[1] += 1
    plt.plot(dates, WWSB      , label = label1)
    plt.plot(dates, WWSB_de   , label = label2)
    plt.xlabel("time")
    plt.ylabel(@sprintf("P [%s]", punit))
    plt.grid()
    plt.legend()
    plt.title("WWSB power detrended")
    plt.savefig(joinpath(fig_dir, path1))

    plt.figure(fig[1]); fig[1] += 1
    plt.plot(dates, ΔEL, label = label3)
    plt.xlabel("time")
    plt.ylabel(@sprintf("(WWSB - LOad) [%s]", punit))
    plt.legend()
    plt.grid()
    plt.title("WWSB power - Load")
    plt.savefig(joinpath(fig_dir, path2))
end

function plot_cumulative_power(dates, WWSB, Load, oprod, punit, fig_dir, fig)
    plt.figure(fig[1]); fig[1] += 1
    cΔEL = cumsum(WWSB .- Load)
    @info 1, cΔEL[end]
    plt.plot(dates, cΔEL, label = "o=1")
    for op in oprod
        cΔEL = cumsum(WWSB.*op .- Load)
        plt.plot(dates, cΔEL, label = @sprintf("o=%f", op))
        @info op, cΔEL[end]
    end
    plt.xlabel("time")
    plt.ylabel(@sprintf("(WWSB - LOad) [%s]", punit))
    plt.legend()
    plt.grid()
    plt.title("cumsum(WWSB power - Load)")
    plt.savefig(joinpath(fig_dir, "cumsumWWSBLoad.png"))
end

function plot_storage_fill_level(dates::Vector{DateTime}, Load_de::Vector{Float64}, WWWSB_de::Vector{Float64}, 
                                 WWSB_scaled::Vector{Vector{Float64}},
                                 storages, oprod, j, fig_dir::String, fig::Vector{Int64}, punit; plot_all_p = false)

    mkpath(fig_dir)
    
    title = @sprintf("S_%d", j)
    pngpath = @sprintf("S_%d.png", j)
    @info title

    prozent = "%"
    plt.figure(fig[1]); fig[1] += 1
    for (stg, op) in zip(storages, oprod)
        label = @sprintf("%s storage_cpacity=%2.fh %s, op = %3.f %s", title, stg.SC, punit, (op-1.0)*100.0, prozent)
        plt.plot(dates, stg.SF, label = label)
    end
    plt.xlabel("time")
    plt.ylabel(@sprintf("storage fill level [%sh]", punit))
    #plt.legend()
    plt.grid()
    plt.savefig(joinpath(fig_dir, pngpath))

    if plot_all_p
        nb_op =length(oprod)
        for i in 1:nb_op
            stg = storages[i]
            op = oprod[i]

            sc = stg.SC

            plt.figure(fig[1]); fig[1] += 1
            plt.plot(dates, WWSB_scaled[i],   label="P")
            plt.plot(dates, Load_de,   label="Load")
            plt.xlabel(@sprintf("WWSB_op [%s]", punit))
            plt.ylabel(@sprintf("Load [%s]", punit))
            plt.title(@sprintf("%s Load, WWSB_op: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            plt.legend()
            plt.savefig(joinpath(fig_dir, @sprintf("%s_1_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            plt.figure(fig[1]); fig[1] += 1
            plt.plot(dates, Load_de,   label="Load")
            plt.plot(dates, stg.I4,   label="IF")
            plt.xlabel(@sprintf("WWSB [%s]", punit))
            plt.ylabel(@sprintf("Load, I4 [%s]", punit))
            plt.title(@sprintf("%s Storage power inflow: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            plt.legend()
            plt.savefig(joinpath(fig_dir, @sprintf("%s_2_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            plt.figure(fig[1]); fig[1] += 1
            plt.plot(dates, Load_de,   label="Load")
            plt.plot(dates, stg.I6,  label="I6")
            plt.xlabel(@sprintf("WWSB [%s]", punit))
            plt.ylabel(@sprintf("Load, I6 [%s]", punit))
            plt.title(@sprintf("%s Storage power outflow: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            plt.legend()
            plt.savefig(joinpath(fig_dir, @sprintf("%s_3_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            plt.figure(fig[1]); fig[1] += 1
            plt.plot(dates, Load_de,   label="Load")
            plt.plot(dates, stg.I7, label="IM")
            plt.xlabel(@sprintf("WWSB [%s]", punit))
            plt.ylabel(@sprintf("Load, I7 [%s]", punit))
            plt.title(@sprintf("%s Power import: SC=%2.f[%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            plt.legend()
            plt.savefig(joinpath(fig_dir, @sprintf("%s_4_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            plt.figure(fig[1]); fig[1] += 1
            plt.plot(dates, Load_de,   label="Load")
            plt.plot(dates, stg.I5, label="CT")
            plt.xlabel(@sprintf("WWSB [%s]", punit))
            plt.ylabel(@sprintf("Load, I5 [%s]", punit))
            plt.title(@sprintf("%s Power curtailment: SC=%2.f %s, op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            plt.legend()
            plt.savefig(joinpath(fig_dir, @sprintf("%s_5_%2.f_%1.f", title, sc, (op-1.0)*100.0)))
        end
    end
end

