include("coefficient_mgr.jl")
# ==============================================================================
# Simple Tests
# ==============================================================================

function SimpleTest()

    cm = CoefficientMatrixManager()
    y = VariableCoefManager(:y, 3,3)
    x = VariableCoefManager(:x, 2,2)
    for idx in 1:3
        y[idx, idx] = 1
    end
    for idx in 1:2
        x[idx, idx] = 1
    end
    RegisterVar(cm, x)
    RegisterVar(cm, y)
    RegisterVarableCoefficients(cm, x)
    RegisterVarableCoefficients(cm, y)
    cm|>GetMatrix|>Matrix|>display
return cm end

SimpleTest()