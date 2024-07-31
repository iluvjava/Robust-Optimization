# Everytime when the experiment is run, this file must be used to setup experiments parameters, 
# and assembles all the files needed to run the main CCGA algorithm. 

using Infiltrator, ProgressMeter, Dates, Distributions, Tables, CSV

include("../src/utilities.jl")

module MatrixConstruct
    "Experiments data folder. "
    DATA_DIR = "data_realistic"

    include("../src/coef_mgr_v2.jl")
    include("../src/problem_parameters.jl")
    
    # Apply multiplers on different parameters for the data. The below is the default profile. 
    CONST_PROBLEM_PARAMETERS.HORIZON = 24
    CONST_PROBLEM_PARAMETERS.Φ = 2400
    for k in keys(DEMAND_RESPONSE.R)
        DEMAND_RESPONSE.R[k] *= 1
    end
    STORAGE_SYSTEM.Capacity[1] = 100*3.5
    PRIMARY_GENERATORS.Pmax .*= 1
    PRIMARY_GENERATORS.Pmin .*= 1
    PRIMARY_GENERATORS.RU .*= 1
    PRIMARY_GENERATORS.RD .*= 1
    @info "Manual multipliers has been applied to the parameters on the data set for experiments. "

    include("../src/matrix_constructions.jl")
    h = rhs

    """
    Change some parameters for the generators ramp-up, down limit, and multipliers for bdges, generation level, 
    demand response. And finally time horizon and budget. The matrices inside of the module will be reconstructed 
    after the parameter change. This should make s  ensitivity analysis easier. 

    """
    function ChangeParametersMakeMatrices(
        time_horizon, 
        budget, 
        demand_response, 
        multiplier, 
        storage_sys_capacity_multiplier, 
        pmax_multiplier, 
        pmin_multiplier, 
        ru_multiplier, 
        rd_multiplier, 
        )
        
        EstablishMatrices()
    end


end


include("../src/ccga_modeling.jl")
include("../src/ccga_loops.jl")



