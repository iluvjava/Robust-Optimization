include("../scratch_paper2/coefficient_mgr.jl")

cmm = CoefficientMatrixManager()
x = VariableCoefManager(:x, 3,3)
y = VariableCoefManager(:y, 2,2)
x[3, 3] = 1
x[1,1] = 1
y[2, 2] = 1

RegisterVar(cmm, x)
RegisterVar(cmm, y)
RegisterVarableCoefficients(cmm, x)
RegisterVarableCoefficients(cmm, y)

cmm |>GetMatrix
