include("../scratch_paper2/coefficient_mgr.jl")
using GLPK, JuMP
model = Model(GLPK.Optimizer)
@variable(model, x[1:5, 1:5], Bin)
@variable(model, y[1:5], Bin)

CoefMatrix1 = rand(5,5)
CoefMatrix2 = rand(5)
x = VariableCoefManager(:x, 5,5)
y = VariableCoefManager(:y, 5)
for II in 1:5, III in 1:5
    x[II, III] = CoefMatrix1[II, III]
end
for II in 1:5
    y[II] = CoefMatrix2[II]
end

