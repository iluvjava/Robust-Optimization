### The point of unit tests is not to catch any bug, but a preliminary verification that we 
### didn't break the code while adding more code to existing code. 
using Test, Infiltrator

include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

"The global environment for the gurobi solver. "
const GUROBI_ENV = Gurobi.Env()

@assert MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON == 4 "Constant TIME_HORIZON is not suitable for unit testing, please change it to: 4. "

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


@testset "Testing Basics" begin
"""
The minimal amount of code to test to make sure the newly introduce code are not breaking the old 
code. 
"""
    
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

    global γ⁺ = 50  # initial parameters
    global d̂ = 100*ones(size(MatrixConstruct.H, 2))
    global w̄    # solved by msp 
    global γ̄    # solved by msp 
    global mp   
    global msp
    global fmp  
    global fsp
    global d⁺   # solved by fmp 
    global q    # solved by fsp

    """
    Verifying that the matrices: C,B,H,G and the vector u, w, q are all having 
    compatible dimensons with each other. 
    """
    function VerifyMatrices()
        @assert w_len == size(B, 2) "vector w and matrix B has incompatible size. "
        @assert u_len == size(C, 2) "Vector u and matrix C has incompatible size. "
        @assert q_len == size(G, 2) "Vector q andmatrix G has incompatible size."
        @assert size(C, 1) == size(B, 1) == size(G, 1) == size(H, 1) "The number of rows for the matrices: C, B, G, H are not equal. "*
        "Currently they are: C=$(size(C, 1)), B=$(size(B, 1)), G=$(size(G, 1)), H=$(size(H, 1)). "
        return true 
    end

    """
    Using the main problem to solve for a configurations to establish the variable 
    w̄, γ̄ for the system. 
    """
    function SetupMainProblem()
        mp = MP(MakeEmptyModel(), γ⁺)
        w̄ = FeasibleConfig(mp, d̂)
        println("solution value w solved by the main problem is: ")
        w|>println
        return w !== nothing
    end
    
    """
    Using the established w̄ solved from the main problem to setup the master problem, several different 
    parameter config will be used to see if the master problem is working ok. 
    """
    function SetupMasterProblem()
        @info "Setting up and solving master problem: "
        msp = MSP(MakeEmptyModel(), d̂, γ⁺, block_demands=0, objective_types=1)
        Solve!(msp)
        γ̄ = msp|>GetGamma
        @info "gamma upper:\n$(γ̄)"
        return !(objective_value(msp)|>isnan)
    end

    """
    Test out the FMP, the FMP with the bilinear reformulations. 
    """
    function SetupFMP()
        @info "Trying to initial solve on the instance of FMP with bilinear reformulations. "
        fmp = FMP(w̄, γ̄, d̂, MakeEmptyModel())
        Solve!(fmp)
        @info "The objective value for initial fmp: $(fmp|>objective_value). "
        d⁺ = GetDemandVertex(fmp)
        return !isnan(fmp|>objective_value)
    end

    """
    Using the constructor of FSP and testing it. 
    """
    function SetupFSP()
        @info "Testing the FSP problem using parameters from the master problem and FMP. "
        fsp = FSP(w̄, d⁺, MakeEmptyModel())
        Solve!(fsp)
        @info "The fsp has objective value of:\n$(fsp|>objective_value). "
        q = Getq(fsp)
        return true
    end
    
    function FSPLessThanFMP()
        return (fsp|>objective_value) <= (fmp|>objective_value)
    end

    function IntroduceCutToFMP()
        @info "Introducing cut to the instance of FMP. "
        IntroduceCut!(fmp, q)
        Solve!(fmp)
        @info "Solving fmp with the cut and the objecive value of fmph is:\n$(fmp|>objective_value). "
        return true
    end

    """
    Using the results of the FMP to add an cut to the master problem and make sure everything is 
    running ok. 
    """
    function IntroduceCutToMSP()
        @info "Introducing the cut back to the master problem and then solving it. "
        IntroduceCut!(msp, GetRhoPlus(fmp), GetRhoMinus(fmp))
        Solve!(msp)
        @info "The gamma upper bound from the msp is: \n$(msp|>objective_value).\nthe gamma vector is:\n $(GetGamma(msp))"
        return true
    end
    
    # actually running these tests. 
    @test VerifyMatrices()
    @test SetupMainProblem()
    @test SetupMasterProblem()
    @test SetupFMP()
    @test SetupFSP()
    @test FSPLessThanFMP()
    @test IntroduceCutToFMP()
    @test IntroduceCutToMSP()


end

w̄ = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, -0.0, 0.0, -0.0, 0.0, 0.0, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, 1.0, -0.0, 0.0, 0.0, -0.0, -0.0, 0.0, -0.0, -0.0, 0.0, 0.0, -0.0, 0.0]
d̂ = 200*ones(size(MatrixConstruct.H, 2))
# γ = [46.06048106344293, 46.06048106344293, 46.06048106344293, 46.06048106344293, 46.06048106344293, 46.06048106344293, 25.528000364915183, 25.528000364915183, 25.528000364915183, 25.528000364915183, 25.528000364915183, 25.528000364915183, 39.570281823383695, 39.570281823383695, 39.570281823383695, 39.570281823383695, 39.570281823383695, 39.570281823383695, 43.71908763816102, 43.71908763816102, 43.71908763816102, 43.71908763816102, 43.71908763816102, 43.71908763816102]
γ = 50*ones(size(MatrixConstruct.H, 2))




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



@testset "FMPH Basic Testing " begin 

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
        return obj_fsp <= obj_value2 + eps(Float32)
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

