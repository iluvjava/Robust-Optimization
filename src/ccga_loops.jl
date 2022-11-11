using Infiltrator, ProgressMeter, Dates
include("utilities.jl")
include("matrix_construction_export.jl")
include("ccga_modeling.jl")


"The global environment for the gurobi solver. "
const GUROBI_ENV = Gurobi.Env()
"The major directories where every session of the ccga results are going to output to. "
const RESULTS_DIRECTORY = "./ccga_results"

if !isdir(RESULTS_DIRECTORY)
    mkdir(RESULTS_DIRECTORY)
end

"""
Get the current system time stamp up to milisec that can be interpolated within files names in most platforms. 
"""
function TimeStamp()
return "["*(Date(now())|>string)*" "*(Time(now())|>string)*"]" end

### SESSION_FILE =======================================================================================================
# stores the location of a file. And it will write to it and close the stream, one line at a time. 
# If writing to the file is frequent then this will be slow. 

"""
Session files models an output files for the CCGA full iterations. 
- `file_loc::String`: The location of the file that we are printing to. 
"""
struct SessionFile
    file_loc::String
    function SessionFile(file_loc::String)
        touch(file_loc) # create
    return new(file_loc) end
end

"""
A () operator for the type SessionFile which **accepts** `String` and opens the file stream then print
the content of the string to the file specified in the field of the type and then 
**return** the string for other processing tasks. 
"""
function (this::SessionFile)(that::String)::String
    open(this.file_loc, "a+") do io
        println(io, that)
    end
return that end

"""
Pass a IO lambada function and apply that lamba to the underlying IOStream of the file. 
"""
function(this::SessionFile)(fxn::Function)::SessionFile
    open(this.file_loc, "a+") do io
        fxn(io)
    end
return this end

"""
Close the IO for the internal file writing stream for the given FileSession instance. 
### Arguments
* `this::SessionFile`: Method is a member method for this type. 
"""
function Close(this::SessionFile)
    close(this.file_stream)
return end


global FILE_SESSTION_TIME_STAMP = replace(replace(TimeStamp(), r":|\."=>"-"), r"\[|\]"=>"")
mkdir(RESULTS_DIRECTORY*"/$FILE_SESSTION_TIME_STAMP")
"The full directry ended without / for storing everything for the CCGA Full Run."
global SESSION_DIR = RESULTS_DIRECTORY*"/"*FILE_SESSTION_TIME_STAMP
"stores the main output for the Full CCGA Routine"
global SESSION_FILE1 = SessionFile(SESSION_DIR*"/"*"main_print_out.txt")
"stores the parameters that are used to run the algorithm."
global SESSION_FILE2 = SessionFile(SESSION_DIR*"/"*"ccga_parameters.txt")
"stores the results of the full CCGA. "
global SESSION_FILE3 = SessionFile(SESSION_DIR*"/"*"ccga_results.txt")         # the results from the ccga outer and inner iterations. 


"""
Make an optimizer. The default is gurobi. It's always gurobi. And all options can be tweaked 
when calling this function which makes the process easier. all the parameters meaning about the 
gurobi solver can be found here [here](https://www.gurobi.com/documentation/9.5/refman/parameters.html). 
# Arguments
- `optimality_gap`: The "MIPGap" for the gurobi solver. 
- `time_out::=1`: Try to terminate the gurobi solver after a certain number of seconds has passed. 
- `solver_name::String`: Give the solver a name, if this parameter is specified then a file with the 
given solver name and a TimeStamp will be stored to the global SESSION_DIR. 
- `mip_focus::Int`: Change the mode of focus for the GUROBI solver. 

"""
function MakeOptimizer(;optimality_gap=0.001, time_out::Int=180, solver_name::String="", mip_focus::Int=0)
    model = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    set_optimizer_attribute(model, "MIPGap", optimality_gap)
    set_optimizer_attribute(model, "TIME_LIMIT", time_out)
    set_optimizer_attribute(model, "MIPFocus", mip_focus)
    if solver_name !== ""
        set_optimizer_attribute(model, "LogFile", SESSION_DIR*"/$(solver_name)_$(TimeStamp())_gurobi_log.txt")
    end
return model end


# ======================================================================================================================
# CCGA Inner Loop struct. 
# ======================================================================================================================

"""
A struct that contains all the parameters involved for the CCGA inner forloop iterations. Let me explain 
all the fields it has: 

### Fields: 
* `fmp`: the last instance of the fmp after the ccga inner forloop exited. 
* `fsp`: the last instance of the fsp after the ccga inner loop exited. 
* `all_ds`: all the d_star ths is being generated by the fmp during the iterations of the inner CCGA loop. 
* `all_qs`: all the q_star provided by the fsp in the CCGA inner loop. 
* `upper_bounds`: The objective values of fmp during the inner interations, stored in a list. 
* `lower_bounds`: The objective values of the fsp during the inner iterations, stored in a list. 
* `terminaton_staus`: A code representing different reasons that made the inner CCGA interminate: 
    * Status code  `0`: terminates because fsp is positive or because fmp and fsp converged to zero. 
    * Status code `-1`: terminates because max iteration counter is reached. 
"""
mutable struct CCGAInnerResults
    fmp::AbsFMP
    fsp::FSP
    all_ds::Vector
    all_qs::Vector
    upper_bounds::Vector
    lower_bounds::Vector
    termination_status::Int 

    function CCGAInnerResults(
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

"""
We report the feasible solutions produced by the FSP. Here is a list of parameters we are 
interested in: 
* `p, p',sr, sr',regu, regu', regd, regd' ,nsp, nsp',g+, g-, rho+, rho-`; 
    And the values of these parameters will be returned as a formatted multi-line text. 

### Arguments
- `this::CCGAInnerResults`: It's a member method of this type. 
"""
function ProduceReport(this::CCGAInnerResults)::String
    string_list = Vector{String}()
    u = this.fsp.u; C = MatrixConstruct.C_
    var_coef_holder_list = [
        MatrixConstruct.p, MatrixConstruct.p′,
        MatrixConstruct.sr, MatrixConstruct.sr′,
        MatrixConstruct.regu, MatrixConstruct.regu′,
        MatrixConstruct.regd, MatrixConstruct.regd′,
        MatrixConstruct.nsp, MatrixConstruct.nsp′,
        MatrixConstruct.g_plus, MatrixConstruct.g_minus, 
    ]
    for var_coef_holder in var_coef_holder_list
        (starting_at, ending_at) = MatrixConstruct.ColumnRegimeFor(C, var_coef_holder)
        for idx = starting_at:ending_at
            push!(string_list, "$(u[idx]) = $(u[idx]|>value)")
        end
    end
    ρ⁺ = reshape(GetRhoPlus(this.fmp), size(MatrixConstruct.d))
    ρ⁻ = reshape(GetRhoMinus(this.fmp), size(MatrixConstruct.d))
    push!(string_list, "ρ⁺:")
    push!(string_list, repr("text/plain", ρ⁺))
    push!(string_list, "ρ⁻: ")
    push!(string_list, repr("text/plain", ρ⁻))
    for idx in length(string_list) + 1 :-1: 1
        insert!(string_list, idx, "\n")
    end
    
return join(string_list) end

"""
Produce a figure that contains the objective values of the FMP and the FSP, plotted on the same graph. 
The x-axis is the number of iterations and the y-axis are the objective values of both FSP, FMP. 
The function will return a figure. 

### Arguments: 
- `this::CCGAInnerResults`: A member method of this type. 
"""
function ProducePlot(this::CCGAInnerResults)::Plots.Plot
    fig = plot(
        (1:length(this.upper_bounds)|>collect) .- 0.5, 
        this.upper_bounds, 
        label="upper_fmp",
        marker=:x
    )
    plot!(
        fig, 
        1:length(this.lower_bounds)|>collect,
        this.lower_bounds,
        marker=:+,
        label="lower_fsp"
    )
return fig end

# ======================================================================================================================
# CCGA Outer Loop struct 
# ======================================================================================================================


"""
### CCGAOuterResults
Another struct that models the data exposed during the iterations of the full CCGA forloop (Inner and Outter). 
It will stores the following items: 

* `inner_loops::Vector{CCGAInnerResults}`: A list of CCGAInnerLoop that is obtained during the iterations of the outer CCGA loops. 
* `msp_objectives::Vector{Float64}`: All the objective values of the msp for each iterations of the CCGA, initial msp without any 
    cut should also be introduced into the vector. 
* `msp_gamma::Vector{Vector{Float64}}`: The vector of gammas from the msp is stored here. 
* `msp::Union{MSP, Nothing}`: final instance of the master problem is stored here, it has all the cuts made for outer CCGA iterations. 
* `termination_status::Int`: 
    - `-1`: max iteration is reached and it didn't break out of the forloop due to any conditions. 
    - `-2`: break out pre-maturally due to inner ccga forloop reaching max iterations counter. 
    - `0`: Successfully finished all iterations without reaching the outer maximum iterations. 
* `fmp_initial_objectives::Vector{Float64}`: The initial objective value for the fmp at the start of each of the inner 
    CCGA for loops. 
* `fsp_initial_objectives::Vector{Floats}`: The initial objective values for the fsp at the start of each of the inner 
    CCGA for loops. 
"""
mutable struct CCGAOuterResults
    inner_loops::Vector{CCGAInnerResults}
    msp_objectives::Vector{Float64}
    msp_gamma::Vector{Vector{Float64}}
    msp::Union{MSP, Nothing}
    termination_status::Int
    fmp_initial_objectives::Vector{Float64}
    fsp_initial_objectives::Vector{Float64}

    function CCGAOuterResults()
        this= new(
            Vector{CCGAInnerResults}(),
            Vector{Float64}(),
            Vector{Vector{Float64}}(),
            nothing,
            0, 
            Vector{Float64}(), 
            Vector{Float64}()
        ) 
    return this end
end

"""
    Add an instance of ccga_inner forloop parameters for the instance. 
"""
function (this::CCGAOuterResults)(that::CCGAInnerResults)::CCGAOuterResults
    push!(this.inner_loops, that)
    push!(this.fmp_initial_objectives, that.upper_bounds[1])
    push!(this.fsp_initial_objectives, that.lower_bounds[1])
return this end

"""
    Call on an instance of `::MSP`, we will records all the gamma of the MSP instance, and then 
    we will store the objective values of the current instance of **MSP**. 
"""
function (this::CCGAOuterResults)(that::MSP)::CCGAOuterResults
    push!(this.msp_gamma, GetGamma(that))
    push!(this.msp_objectives, objective_value(that))
return this end

"""
Returns a *string* that is reporting all the results 
after the finishing the outer forloop iterations. We will be 
reporting the following important parameters: 
1. The initial objective values for the fmp, and fsp during each of the inner CCGA iterations. 
2. what are the upper bound and lower bound from the fmp and fsp during the last CCGA inner iteration. 
3. Make plots and return the figures if it's being asked to do it. 

### Arguments
- `this::CCGAOuterResults`: A member method for this type. 
"""
function ProduceReport(this::CCGAOuterResults)::String
    string_list = Vector{String}()
    push!(string_list, "The msp objective value list is: ")
    push!(string_list, this.msp_objectives|>repr)
    push!(string_list, "The list of gammas it figured out packed into columns of a matrix is: ")
    push!(string_list, repr("text/plain",hcat(this.msp_gamma...)))
    push!(string_list, "The termination status is: $(this.termination_status)")
    for idx in length(string_list) + 1 :-1: 1
        insert!(string_list, idx, "\n")
    end
    
return string_list|>join end


"""
`SaveAllModels(this::CCGAOuterResults)` saves all the the models involved during the execution of the full CCGA 
iterations. The files will be all the final version of *FMP*, *MSP* for the inner and outer iterations in the `MOI`
file format compressed into `.gz`. 
"""
function SaveAllModels(this::CCGAOuterResults)
    for idx in eachindex(this.inner_loops)
        write_to_file(
            this.inner_loops[idx].fmp.model,
            SESSION_DIR*"/"*"fmp_final_inner_$(idx).mof.json.gz",
            )
    end
    write_to_file(this.msp.model, SESSION_DIR*"/"*"msp_final.mof.json.gz")
    return nothing
end


"""
makes plots, there are 2 plots: 
1. All the final objective values of the FSP problems and the 
    final values of the FMP problems. 
2. The FSP and the MSP values for the last inner iterations of the CCGA.  
### Arguments: 
- `this::CCGAOuterResults`: It's a member method of this type. 
"""
function ProducePlots(this::CCGAOuterResults)
    fig1 = this.inner_loops[end]|>ProducePlot
    fig2 = plot(
        this.fmp_initial_objectives, label="fmp_initials_vals"
    )
    plot!(
        fig2, 
        this.fsp_initial_objectives, label="fsp_initial_vals"
    )
return fig1, fig2 end




"""
Performs the CCGA Inter forloops and returns the results for making the cut for the *MSP*. 
## Positional arguments: 
- `gamma_bar::Vector{N1}`: The initial gamma bound for each of the demand decision variable, predicted by the MSP.
- `w_bar::Vector{N2}`: The primary generator setup given by the *MSP*. 
- `d_hat::Vector{N3}`: The center of the demand interval. 
## Named arguments: 
- `epsilon::Float64=0.1`: The absolute tolerance for terminating the inner CCGA for loop. 
- `max_iter::Int=8`: The maximum number of iterations before the for loop terminates. 
## Type Declarations
where {N1<:Number, N2<:Number, N3 <:Number}
"""
function CCGAInnerLoop(
    gamma_bar::Vector{N1}, 
    w_bar::Vector{N2},
    d_hat::Vector{N3};
    epsilon::Float64=0.1,
    max_iter::Int=8
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
    model_fmp = MakeOptimizer(solver_name="FMP", mip_focus=1)
    fmp = FMP(w̄, γ̄, d̂, model_fmp)
    @info "$(TimeStamp()) Inner loop is initialized with fmp, and we are solving the initial fmp. "|>SESSION_FILE1
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
        model_fsp = MakeOptimizer()
        @info "$(TimeStamp()) FSP is made and we are solving it. "|>SESSION_FILE1
        fsp = FSP(w̄, d, model_fsp)
        Solve!(fsp); push!(lowerbound_list, fsp |> objective_value)
        @info "$(TimeStamp()) (FSP Lower, FMP Upper) = $((lowerbound_list[end], upperbound_list[end])) at itr = $II"|>SESSION_FILE1
        if lowerbound_list[end] > ϵ
            @info "$(TimeStamp()) Inner CCGA forloop terminated due to a positive lower of bound of"*
            ": $(lowerbound_list[end]) from FSP which is higher than: ϵ=$(ϵ)."|>SESSION_FILE1
            break
        end
        q = Getq(fsp); push!(all_qs, q)
        Introduce!(fmp, q)
        @info "$(TimeStamp()) CCGA Inner loop continues, constraints is introduced to fmp and we are solving it. "|>SESSION_FILE1
        Solve!(fmp);push!(upperbound_list, objective_value(fmp))
        @assert !(objective_value(fmp)|>isnan) "$(premise) FMP"*
            " is infeasible or unbounded DURING the inner CCGA iterations. "
        if abs(upperbound_list[end] - lowerbound_list[end]) < ϵ
            @info "$(TimeStamp()) Inner CCGA forloop termminated due to convergence of"*
                " FSP and FMP on tolerance level ϵ=$(ϵ), new fmp returns: $(fmp|>objective_value)"|>SESSION_FILE1
            break
        end
        if II == max_iter
            termination_status = -1
        end
    end
    
return CCGAInnerResults(
    fmp,
    fsp,
    all_ds,
    all_qs,
    upperbound_list,
    lowerbound_list, 
    termination_status
) end


"""
Performs the outter forloop of the CCGA algorithm with initialized parameters.

### Positional Arguments
- `d_hat::Vector{N1}`: The center of the demands uncertainty interval. 
- `gamma_upper::N2`: The initial scalar upper bound for all the γ in the uncertainty interval. 

### Keyword Arguments:
- `epsilon`: The tolerance for the lower bound and upper bound between FMP, FSP, and it's used for the termination 
    conditions for the inner forloop. 
- `make_plot`: Whether to make a plot for all the results obtain from the execution of the inner CCGA forloop. 
- `inner_max_itr`: The maximum time s for executing the inner forloop of CCGA. 
- `outer_max_itr`: The maximum number of tiles for excuting the out forloop of CCGA. 
- `make_plot::Bool=true`: Make plots for each of the inner CCGA iterations while running, also make plots after 
    everything is finished. 
- `msp_objective_option::Int=2`: See doc for *MSP* for more. 
- `msp_block_demand_option::Int=1`: See doc for type *MSP* for more information. 
"""
function CCGAOuterLoop(
    d_hat::Vector{N1}, 
    gamma_upper::N2;
    epsilon_inner::N3=0.1, 
    epsilon_outer::N4=0.1,
    inner_max_itr::Int=15,
    outer_max_itr::Int=40,
    make_plot::Bool=true, 
    msp_objective_option::Int=2, 
    msp_block_demand_option::Int=1
) where {N1 <: Number, N2 <: Number, N3 <: Number, N4<:Number}

    context = "During the execution of the outter loop of CCGA: "
    @assert length(d_hat) == size(MatrixConstruct.H, 2) "$(context)"*
    "The length of d_hat is wrong, it's $(length(d_hat)), but we want: $(size(MatrixConstruct.H, 2)). "
    @assert gamma_upper > 0 "$context gamma_upper should be a strictly positive number."
    @assert epsilon_inner > 0 "$context parameter epsilon should be strictly positive but got: $epsilon_inner. "
    @assert inner_max_itr > 0 && outer_max_itr > 0 "$context both the inner_max_itr, outter_max_itr should be larger than zero strictly. " 
    if (epsilon_outer > epsilon_inner)
        @warn "The epsilon tolerance for the outer iteration is strictly less than the inner, this might cause outer forloop iterating indefinitely. "
    end
    # Store the parameters of the CCGA for reports. 
    let
        SESSION_FILE2() do io
            println(io, "d̂: ")
            println(io, repr("text/plain", d_hat))
            println(io, "Gamma upper, or M is: $(gamma_upper)")
            println(io, "epsilon_inner = $epsilon_inner")
            println(io, "epsilon_outer = $epsilon_outer")
            println(io, "inner_max_itr = $inner_max_itr")
            println(io, "outer_max_itr = $outer_max_itr")
            println(io, "HORIZON = $(MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON)")
        end
    end
    # END

    ϵ = epsilon_inner; γ⁺ = gamma_upper; d̂ = d_hat
    model_mp = MakeOptimizer()
    mp = MP(model_mp, γ⁺)
    model_msp = MakeOptimizer(solver_name="MSP")
    msp = MSP(model_msp, d̂, γ⁺, block_demands=msp_block_demand_option, objective_types=msp_objective_option)
    PortOutVariable!(mp, :d) do d fix.(d, d̂, force=true) end
    PortOutVariable!(mp, :v) do v fix.(v, 0, force=true) end
    Solve!(mp)
    w̄ = Getw(mp)
    Solve!(msp)
    γ̄ = GetGamma(msp)

    outer_counter = 0
    outer_results = CCGAOuterResults()
    for _ in 1:outer_max_itr
        outer_counter += 1

        @info "$(TimeStamp()) Outter Forloop itr=$outer_counter" |> SESSION_FILE1
        Results = CCGAInnerLoop(γ̄, w̄, d̂, epsilon=ϵ)
        outer_results(Results)
        
        # SHOULD FACTOR IT OUT AND MAKE IT AS PART OF THE OUTER LOOP STRUCT OBJECTIVE. 
        if make_plot
            fig = Results|>ProducePlot
            fig|>display
        end
        if Results.termination_status == -1
            outer_results.termination_status = -2
            @info "$(TimeStamp()) Outer loop terminated due to inner loop reaching maximum iterations without convergence."|>SESSION_FILE1
            break
        end
        if Results.upper_bounds[end] < epsilon_outer
            @info "$(TimeStamp()) Outer for loop terminated due to convergence of FMP, FSP to an objective value of zero."|>SESSION_FILE1
            break
        end
        @info "$(TimeStamp()) Introduced cut to the msp and we are solving it. "|>SESSION_FILE1
        
        IntroduceCut!(
            msp, 
            GetRhoPlus(Results.fmp), 
            GetRhoMinus(Results.fmp), 
        )

        Solve!(msp)
        outer_results(msp)
        w̄ = Getw(msp)
        γ̄ = GetGamma(msp)
        @info "$(TimeStamp()) Objective value of msp settled at: $(objective_value(msp)). "|>SESSION_FILE1
        @assert !isnan(objective_value(msp)) "$context objective value for msp is NaN. "
        @assert !isinf(objective_value(msp)) "$contex objective value of msp is inf. "
        
    end

    outer_results.msp = msp
    if outer_counter == outer_max_itr 
        outer_results.termination_status = -1
    end

    # Print the results to files! 
    SESSION_FILE3() do io
        write(io, outer_results|>ProduceReport)
        write(io, outer_results.inner_loops[end]|>ProduceReport)
    end
    SaveAllModels(outer_results)
    if make_plot
        fig1, fig2 = outer_results|>ProducePlots
    end
    savefig(fig1, SESSION_DIR*"/"*"last_fmp_fsp")
    savefig(fig2, SESSION_DIR*"/"*"initial_fmp_fsp")
    # END

return outer_results end



### ====================================================================================================================
### Problems we are identifying: 
### Trying to executing the CCGA Inner forloops and see what could be causing the infeasibility to the cut to the master 
### problem. 

ϵ = 0.1
γ_upper = 50
d̂ = 200*(size(MatrixConstruct.H, 2)|>ones)
Results = CCGAOuterLoop(
    d̂,
    γ_upper,
    inner_max_itr=10, 
    outer_max_itr=10, 
    msp_objective_option=2
);

