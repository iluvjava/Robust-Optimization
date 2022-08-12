### We want to know more about the feasibility of the FMP, and we want to know if the corner point demands hypothesis is 
### legit. 

include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

using Debugger, Infiltrator


"""
    A struct that contains all the parameters involved for the CCGA inner forloop iterations. 
"""
mutable struct CCGAInnerloopParameters
    fmp::AbsFMP
    fsp::FSP
    all_ds::Vector
    all_qs::Vector
    upper_bounds::Vector
    lower_bounds::Vector

    function CCGAInnerloopParameters(
        fmp::AbsFMP,
        fsp::FSP,
        all_ds::Vector,
        all_qs::Vector,
        upper_bounds::Vector,
        lower_bounds::Vector
    )
        this = new(fmp, fsp, all_ds, all_qs, upper_bounds, lower_bounds)
    return this end
end


"""
    Performs the CCGA Inter forloops and returns the results for the cut. 
"""
function CCGAInnerLoop(
    gamma_bar::N1, 
    w_bar::Vector{N2}, 
    d_hat::Vector{N3};
    epsilon::Float64=0.1,
    max_iter::Int=8,
    sparse_vee::Bool=true
) where {N1<:Number, N2<:Number, N3 <:Number}
    premise = "during executing CCGA Inner for loop: "
    @assert length(w_bar) == size(MatrixConstruct.B, 2) "$(premis)w̄ has the wrong size. please verify."
    @assert length(d_hat) == size(MatrixConstruct.H, 2) "$(premise)d̂ has the wrong size, please check the code. "
    @assert epsilon >= 0 "$(premise)ϵ for terminating should be non-negative. "
    @assert max_iter >= 0 "$(premise)Maximum iterations for the inner CCGA forloop should be non-negative. "
    @assert gamma_bar >= 0 "$(premise)γ̄ should be non-negative. "
    
    γ̄ = gamma_bar; d̂ = d_hat; ϵ=epsilon; w̄ = w_bar
    lowerbound_list = Vector()
    upperbound_list = Vector()
    all_qs = Vector{Vector}()
    all_ds = Vector{Vector}()
    model_fmp = Model(Gurobi.Optimizer)
    fmp = FMP(w̄, γ̄, d̂, model_fmp, sparse_vee=sparse_vee)
    Solve!(fmp)
    @assert !(objective_value(fmp)|>isnan) "$(premise) FMP is either unbounded or unfeasible on the first solve of FMP. "
    push!(upperbound_list, objective_value(fmp))
    
    fsp = nothing
    for _ in 1:max_iter
        d = GetDemandVertex(fmp); push!(all_ds, d)
        model_fsp = Model(Gurobi.Optimizer)
        fsp = FSP(w̄, γ̄, d, model_fsp, sparse_vee=sparse_vee)
        
        Solve!(fsp); push!(lowerbound_list, fsp |> objective_value)
        q = Getq(fsp); push!(all_qs, q)
        Introduce!(fmp, q)

        Solve!(fmp);push!(upperbound_list, objective_value(fmp))
        
        @assert !(objective_value(fmp)|>isnan) "$(premise) FMP is infeasible or unbounded DURING the inner CCGA iterations. "
        
        if abs(upperbound_list[end] - lowerbound_list[end]) < ϵ
            break
        end
        
    end
    @exfiltrate
    
return CCGAInnerloopParameters(
    fmp,
    fsp,
    all_ds,
    all_qs,
    upperbound_list,
    lowerbound_list
) end

### using the sub routine here. ========================================================================================


