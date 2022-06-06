include("../src/ccga_modeling.jl")
mp = MP()
Solve!(mp)
DebugReport(mp)