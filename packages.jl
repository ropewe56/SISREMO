using Pkg

Pkg.develop(url="/home/wester/Projects/Julia/Packages/BaseUtils.jl")
Pkg.develop(url="/home/wester/Projects/Julia/Packages/PhysConst.jl")

Pkg.add("CSV")
Pkg.add("CurveFit")
Pkg.add("DataFrames")
Pkg.add("DataStructures")
Pkg.add("Dates")
Pkg.add("ForwardDiff")
Pkg.add("HDF5")
Pkg.add("HTTP")
Pkg.add("Interpolations")
Pkg.add("JSON3")
Pkg.add("NLopt")
Pkg.add("OrderedCollections")
Pkg.add("Printf")
Pkg.add("PyPlot")
Pkg.add("Statistics")


################################################

Pkg.rm("BaseUtils")
Pkg.rm("PhysConst")
Pkg.rm("CSV")
Pkg.rm("CurveFit")
Pkg.rm("DataFrames")
Pkg.rm("DataStructures")
Pkg.rm("Dates")
Pkg.rm("ForwardDiff")
Pkg.rm("HDF5")
Pkg.rm("HTTP")
Pkg.rm("Interpolations")
Pkg.rm("JSON3")
Pkg.rm("NLopt")
Pkg.rm("OrderedCollections")
Pkg.rm("Printf")
Pkg.rm("PyPlot")
Pkg.rm("Statistics")
