using Infiltrator, ProgressMeter
include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")


# Global solver settings. 
const GUROBI_ENV = Gurobi.Env()


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
Another struct that models the data exposed during the iterations of the full CCGA forloop. 
"""
mutable struct CCGAResults
    inner_loops::Vector{CCGAInnerloopParameters}
    msp::Union{MSP, Nothing}
    
    function CCGAResults()
    return new(Vector{CCGAInnerloopParameters}(), nothing) end
end

"""
Add an instance of ccga_inner forloop parameters for the instance. 
"""
function (this::CCGAResults)(that::CCGAInnerloopParameters)
    push!(this.inner_loops, that)
return end


"""
Performs the CCGA Inter forloops and returns the results for the cut. 
"""
function CCGAInnerLoop(
    gamma_bar::Vector{N1}, 
    w_bar::Vector{N2}, 
    d_hat::Vector{N3};
    epsilon::Float64=0.1,
    max_iter::Int=8,
    sparse_vee::Bool=true
) where {N1<:Number, N2<:Number, N3 <:Number}
    
    premise = "During executing CCGA Inner for loop: "
    @assert length(w_bar) == size(MatrixConstruct.B, 2) "$(premis)w̄ has the wrong size. please verify."
    @assert length(d_hat) == size(MatrixConstruct.H, 2) "$(premise)d̂ has the wrong size, please check the code. "
    @assert epsilon >= 0 "$(premise)ϵ for terminating should be non-negative. "
    @assert max_iter >= 0 "$(premise)Maximum iterations for the inner CCGA forloop should be non-negative. "
    @assert !(0 in (gamma_bar .>= 0)) "$(premise)γ̄ should be non-negative. "
    
    γ̄ = gamma_bar; d̂ = d_hat; ϵ=epsilon; w̄ = w_bar
    lowerbound_list = Vector{Float64}()
    upperbound_list = Vector{Float64}()
    all_qs = Vector{Vector}()
    all_ds = Vector{Vector}()
    model_fmp = Model(() -> Gurobi.Optimizer(GUROBI_ENV)); set_silent(model_fmp)
    fmp = FMP(w̄, γ̄, d̂, model_fmp, sparse_vee=sparse_vee)
    @info "Inner loop is initialized with fmp, and we are solving the initial fmp. "
    Solve!(fmp)
    if objective_value(fmp)|>isnan
        DebugReport(fmp, "fmp_debug_report_inner_ccga")
        Infiltrator.@exfiltrate
        @assert false "$(premise) FMP is either unbounded or unfeasible on the first solve of FMP. "
    end
    push!(upperbound_list, objective_value(fmp))
    
    fsp = nothing
    for II in 1:max_iter
        d = GetDemandVertex(fmp); push!(all_ds, d)
        model_fsp = Model(() -> Gurobi.Optimizer(GUROBI_ENV)); set_silent(model_fsp)
        fsp = FSP(w̄, d, model_fsp, sparse_vee=sparse_vee)
        Solve!(fsp); push!(lowerbound_list, fsp |> objective_value)
        @info "(FSP Lower, FMP Upper) = $((lowerbound_list[end], upperbound_list[end])) at itr = $II"
        if lowerbound_list[end] > ϵ
            @info "Inner CCGA forloop terminated due to a positive lower of bound of"*
            ": $(lowerbound_list[end]) from FSP which is higher than: ϵ=$(ϵ)."
            break
        end
        q = Getq(fsp); push!(all_qs, q)
        Introduce!(fmp, q)
        @info "CCGA Inner loop continues, constraints is introduced to fmp and we are solving it. "
        Solve!(fmp);push!(upperbound_list, objective_value(fmp))
        @assert !(objective_value(fmp)|>isnan) "$(premise) FMP"*
            " is infeasible or unbounded DURING the inner CCGA iterations. "
        if abs(upperbound_list[end] - lowerbound_list[end]) < ϵ
            @info "Inner CCGA forloop termminated due to convergence of FSP and FMP on tolerance level ϵ=$(ϵ). "
            break
        end
        
    end
    
return CCGAInnerloopParameters(
    fmp,
    fsp,
    all_ds,
    all_qs,
    upperbound_list,
    lowerbound_list
) end


"""
    Performs the outter forloop of the CCGA algorithm with initialized parameters. The list of parameters: 
d̂: 
    the centered which the uncertainty interval is going to be based upon. 
gamma_upper: 
    The initial scalar upper bound for all the γ in the uncertainty interval. 
epsilon: 
    The tolerance for the lower bound and upper bound between FMP, FSP, and it's used for the termination 
    conditions for the inner forloop. 
make_plot: 
    Whether to make a plot for all the results obtain from the execution of the inner CCGA forloop. 
    inner_max_itr: 
        The maximum time s for executing the inner forloop of CCGA. 
    outter_max_itr: 
"""
function CCGAOutterLoop(
    d_hat::Vector{N1}, 
    gamma_upper::N2;
    epsilon::N3=0.1, 
    inner_max_itr::Int=10,
    outter_max_itr::Int=5,
    make_plot::Bool=true
) where {N1 <: Number, N2 <: Number, N3 <: Number}

    context = "During the execution of the outter loop of CCGA: "
    @assert length(d_hat) == size(MatrixConstruct.H, 2) "$(context)"*
    "The length of d_hat is wrong, it's $(length(d_hat)), but we want: $(size(MatrixConstruct.H, 2)). "
    @assert gamma_upper > 0 "$context gamma_upper should be a strictly positive number."
    @assert epsilon > 0 "$context parameter epsilon should be strictly positive but got: $epsilon. "
    @assert inner_max_itr > 0 && outter_max_itr >0 "$context both the inner_max_itr, outter_max_itr should be larger than zero strictly. " 
    ϵ = epsilon; γ⁺ = gamma_upper; d̂ = d_hat

    model_mp = Model(() -> Gurobi.Optimizer(GUROBI_ENV)); mp = MP(model_mp, γ⁺); set_silent(model_mp)
    model_msp = Model(() -> Gurobi.Optimizer(GUROBI_ENV)); msp = MSP(model_msp, d̂, γ⁺); set_silent(model_msp)

    PortOutVariable!(mp, :d) do d fix.(d, d̂, force=true) end
    PortOutVariable!(mp, :v) do v fix.(v, 0, force=true) end
    Solve!(mp)
    w̄ = Getw(mp)
    Solve!(msp)

    OuterCounter = 0
    OuterResults = CCGAResults()
    for II in 1:outter_max_itr
        println("Outter Forloop itr=$II")
        OuterCounter += 1
        Results = CCGAInnerLoop(GetGamma(msp), w̄, d̂, epsilon=ϵ, sparse_vee=false)
        OuterResults(Results)
        if make_plot
            fig = plot(Results.upper_bounds, label="upper_fmp", marker=:x)
            plot!(fig, Results.lower_bounds, label="lower_fsp", marker=:x)
            fig|>display
        end

        if Results.upper_bounds[end] < 1e-4
            @info "Outter for loop terminated due to convergence of FMP FSP to an objective value of zero. "
            break
        end

        IntroduceCut!(
            msp, 
            GetRhoPlus(Results.fmp), 
            GetRhoMinus(Results.fmp)
        )
        @info "Introduced cut to the msp and we are solving it. "
        Solve!(msp)
        @info "Objective value of msp settled at: $(objective_value(msp)). "
        w̄ = Getw(msp)
    end
    OuterResults.msp = msp
return OuterResults end



### ====================================================================================================================
### Problems we are identifying: 
### Trying to executing the CCGA Inner forloops and see what could be causing the infeasibility to the cut to the master 
### problem. 

ϵ = 0.1
γ_upper = 50
d̂ = 200*(size(MatrixConstruct.H, 2)|>ones)
Results = CCGAOutterLoop(d̂, γ_upper);


