using Infiltrator, ProgressMeter, Dates
include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

N = 3*4*6
# fmph1 = FMPH1(zeros(N), 10*ones(24), 100*ones(24))
# Solve!(fmph1)
# fmph2 = FMPH2(zeros(N), 10*ones(24), 100*ones(24), fmph1.lambda[1].|>value)
# Solve!(fmph2)
# d_star = fmph2.d.|>value
fmph_vals = Vector()
fsp_vals = Vector()
d_stars = Vector()
w = [1.0
1.0
1.0
1.0
1.0
1.0
-0.0
-0.0
-0.0
-0.0
0.0
-0.0
-0.0
0.0
0.0
-0.0
-0.0
0.0
0.0
-0.0
0.0
0.0
1.0
0.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
0.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
0.0
0.0
0.0
0.0
0.0
0.0
-0.0
-0.0
-0.0
-0.0
-0.0
-0.0
0.0
0.0
0.0
0.0
1.0
0.0
-0.0
-0.0
0.0
-0.0
-0.0
0.0]

gamma = [
    34.17861212077394
    34.17861212077394
    34.17861212077394
    34.17861212077394
    34.17861212077394
    34.17861212077394
    28.321387879226062
    28.321387879226062
    28.321387879226062
    28.321387879226062
    28.321387879226062
    28.321387879226062
    50.0
    50.0
    50.0
    50.0
    50.0
    50.0
    49.168537538745795
    49.168537538745795
    49.168537538745795
    49.168537538745795
    49.168537538745795
    49.168537538745795]
gamma = round.(gamma)

fmp = FMP(w, gamma, 200*ones(24))
Solve!(fmp)

FmpVal, fmph1, fmph2 = FirstHeuristic!(w, gamma , 200*ones(24), initial_demands=GetDemandVertex(fmp))
push!(fmph_vals, FmpVal)

# d_star = fmph2.d.|>value
# push!(d_stars, d_star)
# fsp = FSP(w, d_star)
# Solve!(fsp)
# push!(fsp_vals, fsp.v|>value)
# q = fsp.q.|>value

# @info "$(TryHeuristic!(fmph1, fmph2, q))"




