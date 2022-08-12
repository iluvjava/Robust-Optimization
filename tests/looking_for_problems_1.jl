### Objectives: 
### 1. Checks why FMP is infeasible for gamma that is overly small. 
### 2. Checks why MSP is infeasible after the cut is introduced. 

include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

"""
    Fix the v slack variables for all constraints that are not related to the system demands. So that the feasibility 
    problem is only testing on whether the demands can be satisfies given some non-negative slacks. 
        * This will be used for the MP, FSP, and FMP. 
        * Mutates the variables for the model. 
"""
function SparsifyVee!(v::Vector{VariableRef}, demand_groups="Demand Balance")
    starting, ending = MatrixConstruct.RHS_Groups[demand_groups]
    for II in setdiff(Set(1:size(MatrixConstruct.H, 1)), Set(starting:ending))
        fix(v[II], 0, force=true)
    end
return nothing end



"""
    Problem: Why does making v sparse causes FMP to be infeasible? 
    Repeatible Experiment 1: 
        1. Set up the default master problem and obtain inigial γ̄, w̄. 
        2. Set up the FMP problem and solve for a specific: d̂, ϵ, M. 
            * Checks if the corresponding FMP instance turns out to be feasible or not. 
    The lower bound for lambda seems to be the problem. Also not sure why there is a lower bound for lambda there. 
"""

let
    ϵ = 0.1
    M = 10
    d̂ = 40*(size(MatrixConstruct.H, 2)|>ones)
    model_msp = Model(
                optimizer_with_attributes(HiGHS.Optimizer, "output_flag" =>true),
            )
    global msp = MSP(
        model_msp, 
        d̂,
        M
    )
    Solve!(msp)
    w̄ = Getw(msp)
    γ̄ = GetGamma(msp)

    model_fmp = Model(Gurobi.Optimizer)
    global fmp = FMP(w̄, γ̄, d̂, model_fmp)
    PortOutVariable!(fmp, :v) do v 
        SparsifyVee!(v[end])
    end

    PortOutVariable!(fmp, :lambda) do λ
        vcat(λ... .|> delete_lower_bound)
    end

    Solve!(fmp)
    if objective_value(fmp) |> isnan
        @warn "FMP is infeasible. γ̄ = $(γ̄)"
    end
end


"""
Problem 2: 
    Why is the master problem infeasible after introducing the cut to it? 
What to do: 
    * Change it to a feasibility check problem with slacks
    * Deleting the constraints and testing when the problem suddenly becomes feasible. 
What to keep: 
    Keep the sparse v, so the conflict is not from all the constant constraints from the cut. 
"""


let

    ϵ = 0.1
    M = 30
    d̂ = 40*(size(MatrixConstruct.H, 2)|>ones)

    model_mp = Model(HiGHS.Optimizer)
    global mp = MP(model_mp, M)
    PortOutVariable!(mp, :d) do d
        fix.(d, d̂, force=true)
    end
    Solve!(mp)
    model_msp = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "output_flag" =>true),
    )
    global msp = MSP(
        model_msp, 
        d̂,
        M
    )
    
    Solve!(msp)
    w̄ = Getw(msp)
    γ̄ = GetGamma(msp)

    model_fmp = Model(Gurobi.Optimizer)
    global fmp = FMP(w̄, γ̄, d̂, model_fmp)
    
    PortOutVariable!(fmp, :v) do v 
        SparsifyVee!(v[end])
    end

    Solve!(fmp)

    U = objective_value(fmp)
    global ρ⁺ = GetRhoPlus(fmp)
    global ρ⁻ = GetRhoMinus(fmp)
    global d = GetDemandVertex(fmp)
    model_fsp =  Model(
        optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
    )
    
    global fsp = FSP(w̄, γ̄, d, model_fsp)
    PortOutVariable!(fsp, :v) do v
        SparsifyVee!(v)
    end
    Solve!(fsp)

    L = objective_value(fsp)
    global u = Getu(fsp)
    global q = Getq(fsp)
    global v = Getv(fsp)
    
    @info "Upper Bound U: $U, Lower Bound L: $L" 
    IntroduceCut!(msp, u, q, ρ⁺, ρ⁻)
    @objective(msp|>GetModel, Min, sum(msp[:s]))


    Solve!(msp)

    if msp |> objective_value |> isnan
        @warn "The master problem is infeasible"
    end

    DebugReport(msp, "msp_after_first_cut")
    DebugReport(mp, "main_problem_referece")
end


### ====================================================================================================================
### Testing FMP feasibility with Dual
### ====================================================================================================================


let 

    ϵ = 0.1
    M = 10
    d̂ = 40*(size(MatrixConstruct.H, 2)|>ones)
    model_msp = Model(
                optimizer_with_attributes(HiGHS.Optimizer, "output_flag" =>true),
            )
    global msp = MSP(
        model_msp, 
        d̂,
        M
    )
    Solve!(msp)
    w̄ = Getw(msp)
    γ̄ = GetGamma(msp)

    model_fmp = Model(Gurobi.Optimizer)
    global fmp = FMP(w̄, γ̄, d̂, model_fmp)
    PortOutVariable!(fmp, :v) do v 
        SparsifyVee!(v[end])
    end
    Solve!(fmp)

end