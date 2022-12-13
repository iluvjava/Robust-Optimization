using Infiltrator, ProgressMeter, Dates
include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

N = 3*4*6

fmph_vals = Vector()
fsp_vals = Vector()
d_stars = Vector()
w = zeros(N)
gamma = 30*ones(24)
d_hat = 200*ones(24)


FmpVal, fmph1, fmph2 = FirstHeuristic!(w, gamma , d_hat)

function AlternatingSolve(fmph1::FMPH1, fmph2::FMPH2)
    objective_vals1 = Vector()
    objective_vals2 = Vector()
    for _ in 1:10
        fmph1.d = fmph2.d.|>value
        Solve!(fmph1)
        push!(objective_vals1, fmph1|>objective_value)
        ChangeLambdas!(fmph2, fmph1.lambda.|>(x)->(x.|>value))
        Solve!(fmph2)
        push!(objective_vals2, fmph2|>objective_value)
    end
    return objective_vals1, objective_vals2
end

vals1 , vals2 = AlternatingSolve(fmph1, fmph2)

# d_star = fmph2.d.|>value
# fsp = FSP(w, d_star)
# Solve!(fsp)
# q = fsp.q.|>value
# eta = TryHeuristic!(fmph1, fmph2, q)

# d_star = fmph2.d.|>value
# fsp = FSP(w, d_star)
# Solve!(fsp)
# q = fsp.q.|>value
# eta = TryHeuristic!(fmph1, fmph2, q)

# d_star = fmph2.d.|>value
# fsp = FSP(w, d_star)
# Solve!(fsp)
# q = fsp.q.|>value
# eta = TryHeuristic!(fmph1, fmph2, q)

# fmph2.dual_con_idx|>println


# ChangeLambdas!(fmph2, fmph2.lambda.|>(x)->(x.|>value))