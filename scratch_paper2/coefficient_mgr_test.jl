include("coefficient_mgr.jl")
# ==============================================================================
# Simple Tests
# ==============================================================================

function SimpleTest()

    cm = CoefficientMatrixManager()
    y = VariableCoefManager(:y, 3,3)
    x = VariableCoefManager(:x, 2,2)
    RegisterVar(cm, x)
    RegisterVar(cm, y)

    # register the coefficients of the matrix for the first constraint
    y[1, 1] = -1
    x[1, 1] = 1
    RegisterVarableCoefficients(cm, x)
    RegisterVarableCoefficients(cm, y)
    
    # register the coefficients of the matrix for the second constraint
    empty!(x); empty!(y)  # must remove previous coefficient or else it accumulates!
    cm|>NextRow   # add a new row for our matrix!
    x[1,1] = 1
    y[2,2] = -1
    
    RegisterVarableCoefficients(cm, x)
    RegisterVarableCoefficients(cm, y)

    # register the coefficients of the matrix for the third constraint
    empty!(x); empty!(y)  # must remove previous coefficient or else it accumulates!
    cm |> NextRow   # add a new row for our matrix!
    x[1,1] = 1
    y[3,3] = -1
    RegisterVarableCoefficients(cm, x)
    RegisterVarableCoefficients(cm, y)

    # register a row that has a lot of coefficients for the variable y. 
    empty!(y)
    cm|>NextRow
    y[1:3, 1:3] .= -2
    RegisterVarableCoefficients(cm ,y)


    # Print out the full matrix.     
    cm|>GetMatrix|>Matrix|>display

return cm end

cm = SimpleTest()