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
    Repeatible Experiment 1: 
        1. Set up the default master problem and obtain inigial γ̄, w̄. 
        2. Set up the FMP problem and solve for a specific: d̂, ϵ, M. 
            * Checks if the corresponding FMP instance turns out to be feasible or not. 
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

    Solve!(fmp)
    if objective_value(fmp) |> isnan
        @warn "FMP is infeasible. γ̄ = $(γ̄)"
    end
end

