module RobustOptim
    include("coef_mgr_v2.jl")
    include("problem_parameters.jl")
    include("./matrix_constructions.jl")
    h = RHS
    export VariableCoefficientHolder
    export CoefficientMatrix
end


