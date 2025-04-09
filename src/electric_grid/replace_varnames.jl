pattern = """prod_CostMWh     Pro_CMWh
C_bato           Bat_CMWho
C_bati           Bat_CMWhi
bat_PowerIn      Bat_Pin
bat_PowerOut     Bat_Pout
bat_Einit        Bat_Einit
bat_ηin          Bat_ηin
bat_ηout         Bat_ηout
C_H2o            H2_CMWho
C_H2i            H2_CMWhi        
H2_PowerIn       H2_Pin
H2_PowerOut      H2_Pout
H2_ηin           H2_ηin       
H2_ηout          H2_ηout      
H2_Einit         H2_Einit     
import_CostMWh   Imp_CMWh
imp_inP          Imp_Pin
export_CostMWh   Exp_CMWh
export_outP      Exp_Pout
curt_CostMWh     Cur_CMWh
res_CostMWh      Res_CMWh"""

pattern2 = """E_load0      E_Load0
E_load       E_Load
E_prod0      E_Pro0
E_prod       E_Pro
E_bato       E_Bato
E_bati       E_Bati
E_H2o        E_H2Oo
E_H2i        E_H2Oi
E_impp       E_Imp
E_res        E_Res
E_expp       E_Exp
E_curt       E_Cur
C_load       C_Load
C_prod       C_Prod
C_impp       C_Imp
C_res        C_Res
C_expp       C_Exp
C_curt       C_Cur"""

function replace_var_names(fname, pattern)
    pa = []
    pattern = split(pattern, "\n")
    for pat in pattern
        pat2 = split(pat)
        push!(pa, (strip(pat2[1]), strip(pat2[2])))
    end

    cont = open(joinpath(@__DIR__, fname*".jl"), "r") do io
        read(io, String)
    end

    for p in pa
        cont = replace(cont, p[1] => p[2])
    end

    cont = open(joinpath(@__DIR__, fname*"_2.jl"), "w") do io
        write(io, cont)
    end
end
split(pattern2, "\n")
replace_var_names("electric_grid_run_2", pattern2)
replace_var_names("electric_grid_2", pattern2)
