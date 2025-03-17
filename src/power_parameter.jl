Base.@kwdef mutable struct PowerParameter
    punit                        = u_GW
    scale_Bio                    = 1.0
    Load_scale                   = 1.0
    Woff_scale                   = 1.0
    Won_scale                    = 1.0
    Solar_scale                  = 1.0
    Bio_scale                    = 1.0
    SF1_factor                   = 1.0
    scale_to_installed_power_p   = true
    scale_with_installed_power_p = false
    averaging_hours              = 24*7
    averaging_method             = [:moving_average, :mean][1]
    plot_p                       = true
    plot_all_p                   = false
    log_p                        = false
end

