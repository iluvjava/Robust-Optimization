# Everytime when the experiment is run, this file must be used to setup experiments parameters, 
# and assembles all the files needed to run the main CCGA algorithm. 

using Infiltrator, ProgressMeter, Dates, Distributions, Tables, CSV

include("../src/utilities.jl")

module MatrixConstruct
    "Experiments data folder. "
    DATA_DIR = "data"

    include("../src/coef_mgr_v2.jl")
    include("../src/problem_parameters.jl")
    
    # Manual changes here, for sensitivity analysis 
    CONST_PROBLEM_PARAMETERS.HORIZON = 12
    CONST_PROBLEM_PARAMETERS.Φ = 6e6
    for k in keys(DEMAND_RESPONSE.R)
        DEMAND_RESPONSE.R[k] *= 1
    end
    STORAGE_SYSTEM.Capacity[1] = 100
    PRIMARY_GENERATORS.Pmax .*= 1
    PRIMARY_GENERATORS.Pmin .*= 1
    PRIMARY_GENERATORS.RU .*= 1
    PRIMARY_GENERATORS.RD .*= 1
    @info "Manual multipliers has been applied to the parameters on the data set for experiments. "

    include("../src/matrix_constructions.jl")
    h = rhs
end


include("../src/ccga_modeling.jl")
include("../src/ccga_loops.jl")



