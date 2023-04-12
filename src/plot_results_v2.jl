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
function plot_powers(dates::Vector{DateTime}, Load_ec::Vector{Float64}, Load::Vector{Float64}, RP_ec::Vector{Float64}, RP::Vector{Float64},
    averaging_hours::Int64, fig_dir::String, punit::String, fig::Vector{Int64})

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, Load_ec, label="Load")
    pl.plot(dates, RP_ec, label="Renewables")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", punit))
    pl.grid()
    pl.title("Energy-Charts data")
    pl.legend()

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, RP, label="Renewables")
    pl.plot(dates, Load, label="Load")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", punit))
    pl.grid()
    pl.title("scaled  data")
    pl.legend()

    if averaging_hours > 1
        pl.title(@sprintf("averaged, %d days window", div(averaging_hours,24)))
        pl.savefig(joinpath(fig_dir, @sprintf("RP_averaged.png")))
    else
        @infoe fig_dir
        pl.savefig(joinpath(fig_dir, @sprintf("RP.png")))
    end
end

function plot_detrended(dates::Vector{DateTime}, RP::Vector{Float64}, RP_de::Vector{Float64}, RP_trend::Vector{Float64}, ΔEL::Vector{Float64},
    Load::Vector{Float64}, Load_de::Vector{Float64}, Load_trend::Vector{Float64},
    punit ::String, fig_dir::String, fig::Vector{Int64}; data_are_averaged = false)

    label1 = "Load"
    label2 = "Load_de"
    path1  = "Load_detrended.png"
    path2  = "Load_diff_detrended.png"
    if data_are_averaged
        label1 = "Load_av"
        label2 = "Load_av_de"
        path1  = "Load_av_detrended.png"
    end

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, Load      , label = label1)
    pl.plot(dates, Load_de   , label = label2)
    pl.plot(dates, Load_trend, "r", linewidth=3, label = "Load_trend")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", punit))
    pl.grid()
    pl.legend()
    pl.title("Load detrended")
    pl.savefig(joinpath(fig_dir, path1))

    label1 = "RP"
    label2 = "RP_de"
    label3 = "RP_de - Load"
    path1  = "RP_detrended.png"
    path2  = "RP_diff_detrended.png"
    if data_are_averaged
        label1 = "RP_av"
        label2 = "RP_av_de"
        label3 = "RP_av_de - Load_av"
        path1  = "RP_av_detrended.png"
        path2  = "RP_av_diff_detrended.png"
    end

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, RP      , label = label1)
    pl.plot(dates, RP_de   , label = label2)
    pl.plot(dates, RP_trend,  "r", linewidth=3, label = "RP_trend")
    pl.xlabel("time")
    pl.ylabel(@sprintf("P [%s]", punit))
    pl.grid()
    pl.legend()
    pl.title("Renewable detrended")
    pl.savefig(joinpath(fig_dir, path1))

    pl.figure(fig[1]); fig[1] += 1
    pl.plot(dates, ΔEL, label = label3)
    pl.xlabel("time")
    pl.ylabel(@sprintf("(P_E - P_L) [%s]", punit))
    pl.legend()
    pl.grid()
    pl.title("Renewable power - Load")
    pl.savefig(joinpath(fig_dir, path2))
end

function plot_storage_fill_level(dates::Vector{DateTime}, L::Vector{Float64}, P::Vector{Float64},
    storages, oprod, j, fig_dir::String, fig::Vector{Int64}, punit ::String; plot_all_p = false)

    title = @sprintf("S_%d", j)
    pngpath = @sprintf("S_%d.png", j)
    @infoe title

    prozent = "%"
    pl.figure(fig[1]); fig[1] += 1
    for (stg, op) in zip(storages, oprod)
        label = @sprintf("%s storage_cpacity=%2.fh %s, op = %3.f %s", title, stg.SC, punit, (op-1.0)*100.0, prozent)
        pl.plot(dates, stg.SF, label = label)
    end
    pl.xlabel("time")
    pl.ylabel(@sprintf("storage fill level [%sh]", punit))
    pl.legend()
    pl.grid()
    pl.savefig(joinpath(fig_dir, pngpath))

    if plot_all_p
        for (stg, op) in zip(storages, oprod)
            sc = stg.SC

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, P .* op ,   label="RP")
            pl.plot(dates, L,   label="Load")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s L, P: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_1_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, L,   label="Load")
            pl.plot(dates, stg.I4,   label="IF")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Storage power inflow: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_2_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, L,   label="Load")
            pl.plot(dates, stg.I6,  label="OF")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Storage power outflow: SC=%2.f [%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_3_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, L,   label="Load")
            pl.plot(dates, stg.I7, label="IM")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Power import: SC=%2.f[%sh], op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_4_%2.f_%1.f", title, sc, (op-1.0)*100.0)))

            pl.figure(fig[1]); fig[1] += 1
            pl.plot(dates, L,   label="Load")
            pl.plot(dates, stg.I5, label="CT")
            pl.xlabel(@sprintf("P [%s]", punit))
            pl.title(@sprintf("%s Power curtailment: SC=%2.f %s, op=%3.f %s", title, sc, punit, (op-1.0)*100.0, prozent))
            pl.legend()
            pl.savefig(joinpath(fig_dir, @sprintf("%s_5_%2.f_%1.f", title, sc, (op-1.0)*100.0)))
        end
    end
end

