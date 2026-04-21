mutable struct Storage
    SC   :: Float64         # storage capacity              [punit_h] TWh
    SF   :: Vector{Float64} # storage fill level            [punit_h] TWh
    ηin  :: Float64         # efficiency while charging
    ηout :: Float64         # efficiency while discharging
    I2   :: Vector{Float64} # direct flow to load           [punit] TW, GW, MW
    I3   :: Vector{Float64} # to storage and/or curtailment [punit] TW, GW, MW
    I4   :: Vector{Float64} # to storage                    [punit] TW, GW, MW
    I5   :: Vector{Float64} # to curtailment                [punit] TW, GW, MW
    I6   :: Vector{Float64} # from storage to load          [punit] TW, GW, MW
    I7   :: Vector{Float64} # other sources to load         [punit] TW, GW, MW
end

function Storage(SC::Float64, nb_steps::Int64; ηin = 1.0, ηout = 1.0)
    SF = zeros(Float64, nb_steps)
    I2 = zeros(Float64, nb_steps)
    I3 = zeros(Float64, nb_steps)
    I4 = zeros(Float64, nb_steps)
    I5 = zeros(Float64, nb_steps)
    I6 = zeros(Float64, nb_steps)
    I7 = zeros(Float64, nb_steps)
    Storage(SC, SF, ηin, ηout, I2, I3, I4, I5, I6, I7)
end


"""
    power_step(stg::Storage, L, P, i)

    L - Load
    P - power
"""
function power_step(stg::Storage, L, P, istep)
    ΔSC = stg.SC - stg.SF[istep-1]
    k = 0
    if P - L > 0.0
        stg.I2[istep] = L
        stg.I3[istep] = P - L
        stg.I4[istep] = min(stg.I3[istep], ΔSC/stg.ηin)
        stg.SF[istep] = stg.SF[istep-1] + stg.I4[istep]*stg.ηin
        stg.I5[istep] = max(0.0, stg.I3[istep]-stg.I4[istep])  # curtailment
        k = 1
    else
        stg.I2[istep] = P
        stg.I6[istep] = min(stg.SF[istep-1] * stg.ηout, L - P)
        stg.SF[istep] = stg.SF[istep-1] - stg.I6[istep]/stg.ηout
        stg.I7[istep] = L - stg.I2[istep] - stg.I6[istep]      # residual load
        k = 2
    end
    k
end
"""
