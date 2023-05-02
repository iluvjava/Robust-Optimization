include("../src/ccga_loops.jl")

global DEMANDS_PROFILES = "data_3/demand_profiles.csv"|>open|>CSV.File
global PROFILE = 20
d̂ = [DEMANDS_PROFILES[PROFILE][idx] for idx in 2:(MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON + 1)]
TOL = 1.0
GAMMA_UPPER = 500

Results = OuterLoop(
    d̂,
    GAMMA_UPPER,
    inner_max_itr=10,
    outer_max_itr=40, 
    objective_types=2,
    inner_epsilon=1, 
    inner_routine=InnerLoopHeuristic, 
    block_demands=0, 
    make_plot=true
);