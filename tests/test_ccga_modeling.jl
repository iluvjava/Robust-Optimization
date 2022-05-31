include("../src/ccga_modeling.jl")

d̂ = 100.0*ones(size(RobustOptim.H, 2))
mp = MP(20)
Solve!(mp)
# CCGA First Iterations: Initializations
w, γ = Getw(mp).|>value, GetGamma(mp).|>value
fmp = FMP(w, γ, d̂)
Solve!(fmp)
d = GetDemandVertex(fmp)
fsp = FSP(w, γ, d)
Solve!(fsp)
q = Getq(fsp)
Introduce!(fmp, q)
Solve!(fmp)
UPPER, LOWER = objective_value(fmp), objective_value(fsp) # This is huge, without a doubt. 

# CCGA Repeated iterations. 



