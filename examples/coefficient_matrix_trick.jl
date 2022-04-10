include("../scratch_paper2/coefficient_mgr.jl")
using GLPK, JuMP
model = Model(GLPK.Optimizer)
@variable(model, x[1:5, 1:5], Bin)
@variable(model, y[1:5], Bin)

CoefMatrix1 = rand(5,5)
CoefMatrix2 = rand(5)
x = VariableCoefManager(:x, 5,5)
y = VariableCoefManager(:y, 5)

# just put in the coefficients for x[i, j] for one constraint, using the 
# Variable CoefManager! 
for II in 1:5, III in 1:5
    x[II, III] = CoefMatrix1[II, III]
end
for II in 1:5
    y[II] = CoefMatrix2[II]
end

cmm = CoefficientMatrixManager()
RegisterVar(cmm, x)
RegisterVar(cmm, y)
RegisterVarableCoefficients(cmm, x)
RegisterVarableCoefficients(cmm, y)

# Sparse Coefficient Matrix!!! 
A = cmm|>GetMatrix

vars = vcat(model[:x][:]..., model[:y][:]...)  # flatten JuMP vars

# AND THEN WE HAVE THE CONSTRAINTS!
@constraint(
    model, c1, A*vars .== 1
)
