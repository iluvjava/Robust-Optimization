include("experiment_context_and_assembling.jl")

# Experiment code here, for example: 
global DEMANDS_PROFILES = "$(MatrixConstruct.DATA_DIR)/demand_profiles.csv"|>open|>CSV.File
global PROFILE = 19
d̂ = [DEMANDS_PROFILES[PROFILE][idx] for idx in 2:(MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON + 1)]
TOL = 1.0
GAMMA_UPPER = 1400

Results = OuterLoop(
    d̂,
    GAMMA_UPPER,
    inner_max_itr=20,
    outer_max_itr=20, 
    objective_types=2,
    inner_epsilon=1.0, 
    inner_routine=InnerLoopHeuristic,
    block_demands=0, 
    make_plot=false, 
    msp_optimality_gap=0.05,
    session_time_out=7200
);

