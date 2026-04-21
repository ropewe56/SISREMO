struct Powers
    Load       :: Vector{Float64}
    Woff       :: Vector{Float64}
    Won        :: Vector{Float64}
    Solar      :: Vector{Float64}
    Bio        :: Vector{Float64}
    WWSB       :: Vector{Float64}

    mean_Load  :: Float64
    mean_Woff  :: Float64
    mean_Won   :: Float64
    mean_Solar :: Float64
    mean_Bio   :: Float64
    mean_WWS   :: Float64

    max_Load   :: Float64 # proxy for installed power
    max_Woff   :: Float64 # proxy for installed power
    max_Won    :: Float64 # proxy for installed power
    max_Solar  :: Float64 # proxy for installed power
    max_Bio    :: Float64 # proxy for installed power
end

function Powers(Load, Woff, Won, Solar, Bio)
    mean_Load = mean(Load)
    mean_Woff = mean(Woff)
    mean_Won  = mean(Won)
    mean_Sol  = mean(Solar)
    mean_Bio  = mean(Bio)
    mean_WWS  = mean_Woff + mean_Won + mean_Sol

    max_Load  = maximum(Load)
    max_Woff  = maximum(Woff)
    max_Won   = maximum(Won )
    max_Solar = maximum(Solar)
    max_Bio   = maximum(Bio )

    WWSB = @. Woff + Won + Solar + Bio

    Powers(Load, Woff, Won, Solar, Bio, WWSB, 
        mean_Load, mean_Woff, mean_Won, mean_Sol, mean_Bio, mean_WWS, 
        max_Load, max_Woff, max_Won, max_Solar, max_Bio)
end

function scale_wind_and_solar_by_load_4(powers::Powers, x)
    WWS = @. powers.Woff * x[2] + powers.Won * x[3] + powers.Solar * x[4]
    scale = (powers.mean_Load - powers.mean_Bio) / mean(WWS) * x[1]
    @. (powers.Woff + powers.Won + powers.Solar) .* scale + Bio
end

function scale_wind_and_solar_by_load_1(powers::Powers, x)
    scale = (powers.mean_Load - powers.mean_Bio) / powers.mean_WWS * x[1]
    @. (power.Woff + power.Won + power.Solar) .* scale + power.Bio
end

function scale_wind_and_solar_bio_by_load_0(powers::Powers, x)
    mean_wwsb = mean(powers.WWSB)
    scale = (powers.mean_Load * x[1]) / mean_wwsb
    powers.WWSB .* scale
end

function get_WWSB_scaled(power::Powers, x, id)
    WWSB = if id == 4
        scale_wind_and_solar_by_load_4(power, x)
    elseif id == 1
        scale_wind_and_solar_by_load_1(power, x)
    elseif id == 0
        scale_wind_and_solar_bio_by_load_0(power, x)
    end
    WWSB
end
