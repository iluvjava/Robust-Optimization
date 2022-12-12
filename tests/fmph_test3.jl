using Infiltrator, ProgressMeter, Dates
include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

N = 3*4*6

fmph_vals = Vector()
fsp_vals = Vector()
d_stars = Vector()
w = zeros(N)
gamma = 10*ones(24)
d_hat = 200*ones(24)


FmpVal, fmph1, fmph2 = FirstHeuristic!(w, gamma , d_hat)

d_star = fmph2.d.|>value
fsp = FSP(w, d_star)
Solve!(fsp)
q = fsp.q.|>value
eta = TryHeuristic!(fmph1, fmph2, q)

d_star = fmph2.d.|>value
fsp = FSP(w, d_star)
Solve!(fsp)
q = fsp.q.|>value
eta = TryHeuristic!(fmph1, fmph2, q)

d_star = fmph2.d.|>value
fsp = FSP(w, d_star)
Solve!(fsp)
q = fsp.q.|>value
eta = TryHeuristic!(fmph1, fmph2, q)

