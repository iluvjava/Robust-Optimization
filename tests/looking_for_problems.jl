### Objectives: 
### 1. Checks why FMP is infeasible for gamma that is overly small. 
### 2. Checks why MSP is infeasible after the cut is introduced. 


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

ϵ = 0.1
M = 40
d̂ = 40*(size(MatrixConstruct.H, 2)|>ones)

