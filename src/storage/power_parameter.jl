root = dirname(dirname(@__DIR__))
get_sisremo_root() = dirname(dirname(@__DIR__))
get_data_dir() = joinpath(root, "data")
get_fig_dir()  = joinpath(root, "figures")

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
    sisremo_root                 = get_sisremo_root()
    data_dir                     = get_data_dir()
    fig_dir                      = get_fig_dir()
    start_year                   = 2016
    end_year                     = 2025
end

