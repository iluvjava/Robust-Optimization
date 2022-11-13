include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

fmpmc = McCormickFMP(true)
Solve!(fmpmc)
DebugReport(fmpmc);