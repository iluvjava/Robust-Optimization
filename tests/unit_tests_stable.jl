### The point of unit tests is not to catch any bug, but a preliminary verification that we 
### didn't break the code while adding more code to existing code. 
using Test

include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

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
        w = FeasibleConfig(mp, d̂)
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
        mp = MSP(MakeEmptyModel(), d̂, γ⁺, block_demands=1, objective_types=1)
        Solve!(mp)
        γ̄ = mp|>GetGamma
        @info "gamma upper:\n$(γ̄)"
        return !(objective_value(mp)|>isnan)
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

    function IntroducingCutToFMP()
        @info "Introducing cut to the instance of FMP. "
        IntroduceCut!(fmp, q)
        Solve!(fmp)
        @info "Solving fmp with the cut"
        return true
    end
    
    # actually running these tests. 
    @test VerifyMatrices()
    @test SetupMainProblem()
    @test SetupMasterProblem()
    @test SetupFMP()
    @test SetupFSP()
    @test FSPLessThanFMP()
    @test IntroducingCutToFMP()


end
