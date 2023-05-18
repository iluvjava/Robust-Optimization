include("../src/ccga_loops.jl")

EXPECTED_DEMANDS = 1800
VARIANCE = 50
GAMMA_UPPER = 500
TOL = 1.0
d̂ = EXPECTED_DEMANDS*(size(MatrixConstruct.H, 2)|>ones)
d̂ += rand(Uniform(-VARIANCE, VARIANCE), d̂|>length)
CSV.write("$SESSION_DIR/d_hat.csv", Tables.table(d̂))
Results = OuterLoop(
    d̂,
    GAMMA_UPPER,
    inner_max_itr=10,
    outer_max_itr=40, 
    objective_types=2,
    inner_epsilon=TOL, 
    inner_routine=InnerLoopMIP, 
    block_demands=0, 
    make_plot=true,
    
);