using Infiltrator, ProgressMeter, Dates
include("utilities.jl")
include("matrix_construction_export.jl")
include("ccga_modeling.jl")


# Global solver settings. 
const GUROBI_ENV = Gurobi.Env()


function TimeStamp()
return "["*(Date(now())|>string)*(Time(now())|>string)*"]" end



# ======================================================================================================================
# CCGA Inner Loop struct. 
# ======================================================================================================================

"""
    A struct that contains all the parameters involved for the CCGA inner forloop iterations. Let me explain 
        all the fields it has: 
    * fmp: the last instance of the fmp after the ccga inner forloop exited. 
    * fsp: the last instance of the fsp after the ccga inner loop exited. 
    * all_ds: all the d_star ths is being generated by the fmp during the iterations of the inner CCGA loop. 
    * all_qs: all the q_star provided by the fsp in the CCGA inner loop. 
    * upper_bounds: The objective values of fmp during the inner interations, stored in a list. 
    * lower_bounds: The objective values of the fsp during the inner iterations, stored in a list. 
    * terminaton_staus: A code representing different reasons that made the inner CCGA interminate: 
        Status code  0: terminates because fsp is positive or because fmp and fsp converged to zero. 
        Status code -1: terminates because max iteration counter is reached. 
"""
mutable struct CCGAInnerloop
    fmp::AbsFMP
    fsp::FSP
    all_ds::Vector
    all_qs::Vector
    upper_bounds::Vector
    lower_bounds::Vector

    termination_status::Int 

    function CCGAInnerloop(
        fmp::AbsFMP,
        fsp::FSP,
        all_ds::Vector,
        all_qs::Vector,
        upper_bounds::Vector,
        lower_bounds::Vector, 
        termination_status
    )
        this = new(fmp, fsp, all_ds, all_qs, upper_bounds, lower_bounds, termination_status)
    return this end
end

# ======================================================================================================================
# CCGA Outer Loop struct 
# ======================================================================================================================


"""
    Another struct that models the data exposed during the iterations of the full CCGA forloop. It will stores the 
    following items: 
    inner_loops: A list of CCGAInnerLoop that is obtained during the iterations of the outer CCGA loops. 
    msp_objectives: All the objective values of the msp for each iterations of the CCGA, initial msp without any 
        cut should also be introduced into the vector. 
    msp_gamma: The vector of gammas from the msp is stored here. 
    msp: final instance of the master problem is stored here, it has all the cuts made for outer CCGA iterations. 
    termination_status: 
        status code -1: max iteration is reached and it didn't break out of the forloop due to any conditions. 
        status code -2: break out pre-maturally due to inner ccga forloop reaching max iterations counter. 
        status code  0: Successfully finished all iterations without reaching the outer maximum iterations. 
        
"""
mutable struct CCGAOuterLoop
    inner_loops::Vector{CCGAInnerloop}
    msp_objectives::Vector{Float64}
    msp_gamma::Vector{Vector{Float64}}
    msp::Union{MSP, Nothing}
    termination_status::Int
    
    function CCGAOuterLoop()
    return new(
        Vector{CCGAInnerloop}(),
        Vector{Float64}(),
        Vector{Vector{Float64}}(),
        nothing,
        0
    ) end
end

"""
    Add an instance of ccga_inner forloop parameters for the instance. 
"""
function (this::CCGAOuterLoop)(that::CCGAInnerloop)
    push!(this.inner_loops, that)
return this end

"""
    Given an instance of the msp, we record 
"""
function (this::CCGAOuterLoop)(that::MSP)
    push!(this.msp_gamma, GetGamma(that))
    push!(this.msp_objectives, objective_value(that))
return this end

"""
    return a string that is reporting all the results after the finishing the outer forloop iterations. We will be 
    reporting the following important parameters: 
"""
function ProduceReport(this::CCGAOuterLoop)

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
    @info "$(TimeStamp()) Inner loop is initialized with fmp, and we are solving the initial fmp. "
    Solve!(fmp)
    if objective_value(fmp)|>isnan
        DebugReport(fmp, "fmp_debug_report_inner_ccga")
        @assert false "$(premise) FMP is either unbounded or unfeasible on the first solve of FMP. "
    end
    push!(upperbound_list, objective_value(fmp))
    
    fsp = nothing
    termination_status = 0
    for II in 1:max_iter
        d = GetDemandVertex(fmp); push!(all_ds, d)
        model_fsp = Model(() -> Gurobi.Optimizer(GUROBI_ENV)); set_silent(model_fsp)
        @info "$(TimeStamp()) FSP is made and we are solving it. "
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
        @info "$(TimeStamp()) CCGA Inner loop continues, constraints is introduced to fmp and we are solving it. "
        Solve!(fmp);push!(upperbound_list, objective_value(fmp))
        @assert !(objective_value(fmp)|>isnan) "$(premise) FMP"*
            " is infeasible or unbounded DURING the inner CCGA iterations. "
        if abs(upperbound_list[end] - lowerbound_list[end]) < ϵ
            @info "Inner CCGA forloop termminated due to convergence of FSP and FMP on tolerance level ϵ=$(ϵ). "
            break
        end
        if II == max_iter
            termination_status = -1
        end
    end
    
return CCGAInnerloop(
    fmp,
    fsp,
    all_ds,
    all_qs,
    upperbound_list,
    lowerbound_list, 
    termination_status
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
    epsilon_inner::N3=0.1, 
    inner_max_itr::Int=10,
    outer_max_itr::Int=10,
    make_plot::Bool=true, 
    smart_cut::Bool=true
) where {N1 <: Number, N2 <: Number, N3 <: Number}

    context = "During the execution of the outter loop of CCGA: "
    @assert length(d_hat) == size(MatrixConstruct.H, 2) "$(context)"*
    "The length of d_hat is wrong, it's $(length(d_hat)), but we want: $(size(MatrixConstruct.H, 2)). "
    @assert gamma_upper > 0 "$context gamma_upper should be a strictly positive number."
    @assert epsilon_inner > 0 "$context parameter epsilon should be strictly positive but got: $epsilon_inner. "
    @assert inner_max_itr > 0 && outer_max_itr >0 "$context both the inner_max_itr, outter_max_itr should be larger than zero strictly. " 
    ϵ = epsilon_inner; γ⁺ = gamma_upper; d̂ = d_hat

    model_mp = Model(() -> Gurobi.Optimizer(GUROBI_ENV)); set_silent(model_mp)
    mp = MP(model_mp, γ⁺)
    model_msp = Model(() -> Gurobi.Optimizer(GUROBI_ENV)); set_silent(model_msp)
    msp = MSP(model_msp, d̂, γ⁺)

    PortOutVariable!(mp, :d) do d fix.(d, d̂, force=true) end
    PortOutVariable!(mp, :v) do v fix.(v, 0, force=true) end
    Solve!(mp)
    w̄ = Getw(mp)
    Solve!(msp)
    γ̄ = GetGamma(msp)
    Σγ = objective_value(msp)

    OuterCounter = 0
    OuterResults = CCGAOuterLoop()
    for _ in 1:outer_max_itr
        OuterCounter += 1
        println("Outter Forloop itr=$OuterCounter")
        Results = CCGAInnerLoop(γ̄, w̄, d̂, epsilon=ϵ, sparse_vee=false)
        OuterResults(Results)
        
        # SHOULD FACTOR IT OUT AND MAKE IT AS PART OF THE OUTER LOOP STRUCT OBJECTIVE. 
        if make_plot
            fig = plot(Results.upper_bounds, label="upper_fmp", marker=:x)
            plot!(fig, Results.lower_bounds, label="lower_fsp", marker=:x)
            fig|>display
        end

        if Results.termination_status == -1
            OuterResults.termination_status = -2
            @info "Outer loop terminated due to inner loop reaching maximum iterations without convergence."
            break
        end

        if Results.upper_bounds[end] < 1e-4
            @info "Outer for loop terminated due to convergence of FMP, FSP to an objective value of zero."
            break
        end

        @info "$(TimeStamp()) Introduced cut to the msp and we are solving it. "
        
        IntroduceCut!(
            msp, 
            GetRhoPlus(Results.fmp), 
            GetRhoMinus(Results.fmp), 
            Σγ
        )
        Solve!(msp)
        OuterResults(msp)  
        w̄ = Getw(msp)
        γ̄ = GetGamma(msp)
        @info "Objective value of msp settled at: $(objective_value(msp)). "
        @assert !isnan(objective_value(msp)) "$context objective value for msp is NaN. "
        @assert !isinf(objective_value(msp)) "$contex objective value of msp is inf. "
        if objective_value(msp) < Σγ && OuterCounter > 1 && smart_cut
            Σγ = objective_value(msp)
            DeleteAllPreviousCut!(msp)
            @info "SmartCut is deleting all previous cut due to strict decrease of the msp objective. "
        end
        
    end

    OuterResults.msp = msp
    if OuterCounter == outer_max_itr 
        OuterResults.termination_status = -1
    end
return OuterResults end



### ====================================================================================================================
### Problems we are identifying: 
### Trying to executing the CCGA Inner forloops and see what could be causing the infeasibility to the cut to the master 
### problem. 

ϵ = 0.1
γ_upper = 50
d̂ = 200*(size(MatrixConstruct.H, 2)|>ones)
Results = CCGAOutterLoop(d̂, γ_upper, smart_cut=true, inner_max_itr=10, outer_max_itr=20);


