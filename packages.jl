using Pkg

Pkg.develop(url="/home/wester/Projects/Julia/Packages/BaseUtils")
Pkg.develop(url="/home/wester/Projects/Julia/Packages/PhysConst")

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
Pkg.add("Printf")
Pkg.add("Statistics")

