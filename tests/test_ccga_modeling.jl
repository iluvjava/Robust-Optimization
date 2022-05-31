include("../src/ccga_modeling.jl")

mp = MP()
Solve!(mp)
w = Getw(mp).|>value
γ = GetGamma(mp).|>value
