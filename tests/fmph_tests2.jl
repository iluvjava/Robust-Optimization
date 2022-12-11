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
w = zeros(N)
gamma = zeros(24)
d_hat = 200*ones(24)

FmpVal, fmph1, fmph2 = FirstHeuristic!(w, gamma , d_hat)
fmp = FMP(w, gamma, d_hat)
Solve!(fmp)
