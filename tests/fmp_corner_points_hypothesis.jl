### We want to know more about the feasibility of the FMP, and we want to know if the corner point demands hypothesis is 
### legit. 
include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")


### Testing the FMP directly with Inner CCGA iterations. 

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
M = 16
d̂ = 100*(size(MatrixConstruct.H, 2)|>ones)
global lowerbound_list = Vector()
global upperbound_list = Vector()
global all_qs = Vector{Vector}()
global all_ds = Vector{Vector}()
global d_continuous = Vector{Vector}()

model_mp = Model(HiGHS.Optimizer); global mp = MP(model_mp, M)
PortOutVariable!(mp, :d) do d
    fix.(d, d̂, force=true)
end
PortOutVariable!(mp, :v) do v
    fix.(v, 0, force=true)
end
Solve!(mp)
global w̄ = Getw(mp)
γ̄ = M
model_fmp = Model(Gurobi.Optimizer)
global fmp = FMP(w̄, γ̄, d̂, model_fmp)
Solve!(fmp)
@assert !(objective_value(fmp)|>isnan)

# push!(upperbound_list, objective_value(fmp))
for III in 1:10
    global d = GetDemandVertex(fmp); push!(all_ds, d); push!(d_continuous, fmp.d.|>value)

    model_fsp = Model(Gurobi.Optimizer)
    global fsp = FSP(w̄, γ̄, d, model_fsp)
    
    # PortOutVariable!(fsp, :v) do v
    #     SparsifyVee!(v)
    # end
    
    Solve!(fsp)
    global q = Getq(fsp); push!(all_qs, q)
    Introduce!(fmp, q)
    Solve!(fmp)
    push!(lowerbound_list, fsp |> objective_value)
    @assert !(objective_value(fmp)|>isnan)
    push!(upperbound_list, objective_value(fmp))
    if abs(upperbound_list[end] - lowerbound_list[end]) < 0.1
        break
    end
end
fig = plot(upperbound_list, label="upper_fmp")
plot!(fig, lowerbound_list, label="lower_fsp")
fig|>display

msp_model = Model(Gurobi.Optimizer)
msp = MSP(msp_model, d̂, M)
Solve!(msp)
IntroduceCut!(msp, Getu(fsp), Getq(fsp), GetRhoPlus(fmp), GetRhoMinus(fmp))
Solve!(msp)
DebugReport(msp)
# yeah definitely not that feasible. 

