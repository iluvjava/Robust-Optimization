include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")
using Test, Infiltrator

w̄ = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, -0.0, 0.0, -0.0, 0.0, 0.0, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, 1.0, -0.0, 0.0, 0.0, -0.0, -0.0, 0.0, -0.0, -0.0, 0.0, 0.0, -0.0, 0.0]
d̂ = 200*ones(size(MatrixConstruct.H, 2))
# γ = [46.06048106344293, 46.06048106344293, 46.06048106344293, 46.06048106344293, 46.06048106344293, 46.06048106344293, 25.528000364915183, 25.528000364915183, 25.528000364915183, 25.528000364915183, 25.528000364915183, 25.528000364915183, 39.570281823383695, 39.570281823383695, 39.570281823383695, 39.570281823383695, 39.570281823383695, 39.570281823383695, 43.71908763816102, 43.71908763816102, 43.71908763816102, 43.71908763816102, 43.71908763816102, 43.71908763816102]
γ = 50*ones(size(MatrixConstruct.H, 2))


"The global environment for the gurobi solver. "
const GUROBI_ENV = Gurobi.Env()

"""
    MakeOptimizer(;optimality_gap=0.001, time_out::Int=180, solver_name::String="", mip_focus::Int=0)

Make an gurobi optimizer with all the correct default settings. 
"""
function MakeEmptyModel(;optimality_gap=0.001, time_out::Int=180, solver_name::String="", mip_focus::Int=0)
    model = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    set_optimizer_attribute(model, "LogToConsole", 0)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "MIPGap", optimality_gap)
    set_optimizer_attribute(model, "TIME_LIMIT", time_out)
    set_optimizer_attribute(model, "MIPFocus", mip_focus)
    
    if solver_name !== ""
        set_optimizer_attribute(model, "LogFile", SESSION_DIR*"/$(solver_name)_$(TimeStampConvert)_gurobi_log.txt")
    end
return model end


C = MatrixConstruct.C
B = MatrixConstruct.B
H = MatrixConstruct.H
G = MatrixConstruct.G
h = MatrixConstruct.h
u = MatrixConstruct.u
u_len = u.|>length|>sum
w = MatrixConstruct.w
w_len = w.|>length|>sum
q = MatrixConstruct.q
q_len = q.|>length|>sum


@testset "FMPH Basic Testing " begin 
    
    # testing parameters, for testing things in the interactive REPL. 
    global fmph1
    global fmph2
    global fmphs
    global fsp
    global d_star
    
    function FMPHBasicRun()
        @info "Setting up the FMPH1, 2 with some basic parameters. We are looking at FSP, and the objective value of FMPH2. "
        _, obj_value1, fmph1, fmph2 = FirstHeuristic!(w̄, γ, d̂, MakeEmptyModel)
        d1 = GetDemands(fmph2)
        l1 = GetLambdas(fmph1)
        fmph1 = RebuildFMPH1(fmph1, GetDemands(fmph2), model=MakeEmptyModel())
        _, obj_value2 = TryHeuristic!(fmph1, fmph2)
        d2 = GetDemands(fmph2)
        l2 = GetLambdas(fmph1)
        fmp = FMP(w̄, γ, d̂, MakeEmptyModel())
        Solve!(fmp)
        obj_value3 = fmp|>objective_value
        @info "2 runs of objective value: [$obj_value1, $obj_value2]"
        @info "Classic discrete FMP objective is $obj_value3"
        "2 runs of the lambdas are the same: $((abs.(hcat(l1...) - hcat(l2...)))|>sum <= Float64|>eps)"|>println
        "2 runs of the demands are the same: $((abs.(d1 - d2))|>sum <= Float64|>eps)"|>println

        d_star = GetDemands(fmph2)
        fsp = FSP(w̄, d_star, MakeEmptyModel())
        Solve!(fsp)
        @info "q value obtained for FSP. "
        Getq(fsp)|>println
        @info "objective of fsp: "
        obj_fsp = fsp|>objective_value
        obj_fsp|>println
        @info "objective of fmph2 heuristic"
        obj_value2|>println
        return obj_fsp <= obj_value2
    end

    """
    Test the basic functionality of the FMPHstepper. 
    """
    function FMPHStepperBasic()
        """
        Testing out the functionality of the stepper for the FMPH instances. 
        """|>println
        fmphs = FMPHStepper(w̄, γ, d̂, MakeEmptyModel)
        fmphs()
        fmphs()
        local fsp = FSP(w̄, GetDemands(fmphs), MakeEmptyModel())
        Solve!(fsp)
        q = Getq(fsp)
        fmphs(q)
        fmphs()
        @info "The list of demands for the FMPH stepper is: $(fmphs.obj2). \nThe objective value of fsp was: $(fsp|>objective_value). "
        return true
    end

    """
    1. perform a bilinear hearustic search. 
    2. Change the demand vector and then perform that again. 
    3. Add a vector q from the FSP to it and then perform bilinear search. 
    4. Change the value of demand vector and then perform it again. 
    """
    function FMPHStepperFull()
        fmphs = FMPHStepper(w̄, γ, d̂, MakeEmptyModel)
        fmphs()
        fmphs()
        TryNewDemand(fmphs)
        fmphs()
        fmphs()
        local fsp = FSP(w̄, GetDemands(fmphs), MakeEmptyModel())
        Solve!(fsp)
        q = Getq(fsp)
        fmphs(q)
        TryNewDemand(fmphs)
        TryNewDemand(fmphs)
        TryNewDemand(fmphs)
        @info "This is the list of fmph2 objectis for the fmphs full test:\n$(fmphs.obj2)"
        return true
    end

    @test FMPHBasicRun()
    @test FMPHStepperBasic()
    @test FMPHStepperFull()

end 
