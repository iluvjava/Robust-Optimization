include("experiment_context_and_assembling.jl")

# Experiment code here, for example: 
global DEMANDS_PROFILES = "$(MatrixConstruct.DATA_DIR)/demand_profiles.csv"|>open|>CSV.File
global PROFILE = 42
d̂ = [DEMANDS_PROFILES[PROFILE][idx] for idx in 2:(MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON + 1)]
TOL = 1.0
GAMMA_UPPER = 1000

Results = OuterLoop(
    d̂,
    GAMMA_UPPER,
    inner_max_itr=50,
    outer_max_itr=20,
    objective_types=2,
    epsilon=0.1, 
    inner_routine=InnerLoopHeuristic, 
    # inner_routine=InnerLoopMIP, 
    block_demands=0, 
    make_plot=false, 
    msp_optimality_gap=0.01, 
    session_time_limit=36000
);