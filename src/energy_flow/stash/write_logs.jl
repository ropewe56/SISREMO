function write_to_log(k, istep, j, load, scaled_renewables, mode)
    out1 = open(joinpath(@__DIR__, "log1.log"), mode)
    out2 = open(joinpath(@__DIR__, "log2.log"), mode)

    write(out1, @sprintf("%3d, %d, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n", 
                istep, k, load, scaled_renewables,
                storages[j].I2[i],
                storages[j].I3[i],
                storages[j].I4[i],
                storages[j].I5[i],
                storages[j].I6[i],
                storages[j].I7[i],
                storages[j].SF[i],
                ))
    write(out2, @sprintf("%3d, %d, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n", 
                istep, k2, storages[j-1].I7[istep], storages[j-1].I5[istep],
                storages[j].I2[istep],
                storages[j].I3[istep],
                storages[j].I4[istep],
                storages[j].I5[istep],
                storages[j].I6[istep],
                storages[j].I7[istep],
                storages[j].SF[istep],
                ))
end

function write_to_log(stg::Storage, ΔPin, out, i, j)
    B = stg.IF[i] - stg.I6[i] - (stg.SF[i] - stg.SF[i-1])
    write(out, @sprintf("%5d  %d  ΔP = %9.2e, S = %9.2e, I = %9.2e, O = %9.2e, C = %9.2e, M = %9.2e, B = %9.2e\n", i, j, ΔPin, stg.SF[i], stg.IF[i], stg.OF[i], stg.CT[i], stg.RE[i], B))
end
function write_power_step_to_log(stg::Storage, out, i, j)
    write(out, @sprintf("%5d  %d  P = %9.2e, L = %9.2e, I2 = %9.2e, I3 = %9.2e, I4 = %9.2e, I5 = %9.2e, I6 = %9.2e, I7 = %9.2e\n",
        i, j, stg.SF[i], P, L, stg.I2[i], stg.I3[i], stg.I4[i], stg.I5[i], stg.I6[i], stg.I7[i]))
end
