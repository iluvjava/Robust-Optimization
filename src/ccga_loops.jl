"The global environment for the gurobi solver. "
const GUROBI_ENV = Gurobi.Env()
"The major directories where every session of the ccga results are going to output to. "
const RESULTS_DIRECTORY = "./ccga_results"
"The parameters for default time out options for the solver made for JuMP Models. "
const SOLVER_TIME_OUT = 1800
"The total amount of time allowed for executing the CCGA algorithm. "
global SESSION_TIME_OUT = 1800

if SESSION_TIME_OUT <= 0
    error("SESSION_TIME_OUT out can't be <= 0. ")
end

"The int epoch time for when the session was started. This can be modified. "
global SESSION_START_TIME = nothing


if !isdir(RESULTS_DIRECTORY)
    mkdir(RESULTS_DIRECTORY)
end

"""
function returns the number of seconds since `SESSION_START_TIME`. If -1 
is returned, that means `SESSION_START_TIME` was never established before 
this function is being called. 

"""
function Elapsed()::Int
    if SESSION_START_TIME === nothing
        return -1
    end
    return (now()|>datetime2unix|>floor|>Int) - SESSION_START_TIME
end

"""
Change the solver timeout to be the amount of time remains for this section. 
"""
function AdaptSolverTimeout(model::Model)::Model
    if Elapsed() == -1
        return model
    end
    timeRemains = SESSION_TIME_OUT - Elapsed()
    if timeRemains <= 0 
        error("We run out of time, forced terminations by errors. ")
    end
    set_optimizer_attribute(model, "TIME_LIMIT", timeRemains)
    return model
end

"""
Get the current system time stamp up to mili sec that can be interpolated within files names in most platforms. 
"""
function TimeStamp()
return "["*(Date(now())|>string)*" "*(Time(now())|>string)*"]" end

"""
Convert the timestamp of current time to a format that is file name friendly. 
"""
function TimeStampConvert()
    return replace(replace(TimeStamp(), r":|\."=>"-"), r"\[|\]"=>"")
end

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
function (this::SessionFile)(fxn::Function)::SessionFile
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


# Session file handers are global variables ============================================================================

global FILE_SESSTION_TIME_STAMP = TimeStampConvert()
mkdir(RESULTS_DIRECTORY*"/$FILE_SESSTION_TIME_STAMP")

"The full directry ended without / for storing everything for the CCGA Full Run."
global SESSION_DIR = RESULTS_DIRECTORY*"/"*FILE_SESSTION_TIME_STAMP
"stores the main output for the Full CCGA Routine"
global SESSION_FILE1 = SessionFile(SESSION_DIR*"/"*"main_print_out.txt")
"stores the parameters that are used to run the algorithm."
global SESSION_FILE2 = SessionFile(SESSION_DIR*"/"*"ccga_parameters.txt")
"It stores results at the end of CCGA, which are all the decision variables from the last instance of FSP. "
global SESSION_FILE3 = SessionFile(SESSION_DIR*"/"*"ccga_results.txt")         # the results from the ccga outer and inner iterations. 


"""
Make an instance of Empty JuMP model with pre-attatched optimizers. 
The default is gurobi. It's always gurobi. And all options can be tweaked 
when calling this function which makes the process easier. all the parameters meaning about the 
gurobi solver can be found here [here](https://www.gurobi.com/documentation/9.5/refman/parameters.html). 

# Named Arguments
- `optimality_gap`: The "MIPGap" for the gurobi solver. 
- `time_out::=1`: Try to terminate the gurobi solver after a certain number of seconds has passed. 
- `solver_name::String`: Give the solver a name, if this parameter is specified then a file with the 
given solver name and a TimeStamp will be stored to the global SESSION_DIR. 
- `mip_focus::Int`: Change the mode of focus for the GUROBI solver. 
- `log_to_console::Int=0`: Whether to log the JuMP printout to the console. 
"""
function GetJuMPModel(
    ;optimality_gap=0.05, 
    time_out::Int=SOLVER_TIME_OUT, 
    solver_name::String="", 
    mip_focus::Int=0, 
    log_to_console::Int=0
)
    model = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    set_optimizer_attribute(model, "MIPGap", optimality_gap)
    set_optimizer_attribute(model, "TIME_LIMIT", time_out)
    set_optimizer_attribute(model, "MIPFocus", mip_focus)
    set_optimizer_attribute(model, "LogToConsole", log_to_console)
    if solver_name !== ""
        set_optimizer_attribute(model, "LogFile", SESSION_DIR*"/$(solver_name)_$(TimeStampConvert)_gurobi_log.txt")
    end
return model end



# ======================================================================================================================
# CCGA Inner Loop struct. 
# ======================================================================================================================

"""
CCGAIR: CCGA Innerloop Results
A type that should inherit the following methods:
- `ProduceReport(this::IRBHeuristic)::String`
- `ProducePlot(this::IRBHeuristic)::Plots.Plot`
"""
abstract type CCGAIR

end



"""
CCGAIRBR: CCGA Inner Results Binlinear Reformulations (via MIP).

    CCGAIRBR(
        fmp::AbsFMP,
        fsp::FSP,
        all_ds::Vector,
        all_qs::Vector,
        upper_bounds::Vector,
        lower_bounds::Vector, 
        termination_status
    )

A struct that contains all the parameters involved for the CCGA inner forloop iterations. Let me explain 
all the fields it has: 

### Fields: 
* `fmp`,`fsp`, `all_ds`, `all_qs`, `upper_bounds`, `lower_bounds`, `terminaton_staus`, 
"""
mutable struct IRBReform <: CCGAIR
    
    "the last instance of the fmp after the ccga inner forloop exited. "
    fmp::AbsFMP
    "the last instance of the fsp after the ccga inner loop exited. "
    fsp::FSP
    "all the d_star ths is being generated by the fmp during the last iterations of the inner CCGA loop. "
    all_ds::Vector
    "all the q_star provided by the fsp in the CCGA inner loop. "
    all_qs::Vector
    "The objective values of fmp during the inner interations, stored in a list. "
    upper_bounds::Vector
    "The objective values of the fsp during the inner iterations, stored in a list. "
    lower_bounds::Vector
    """
    A code representing different reasons that made the inner CCGA interminate: 
    * Status code  `1`: terminates because fsp is positive, demand is good for introducing cut on MSP. 
    * status code  `0`: fmp and fsp both converged to zero. 
    * Status code `-1`: terminates because max iteration counter is reached. 
    * Status code `-2`: terminates due to session time out. 
    """
    termination_status::Int 

    function IRBReform(
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
CCGAIRBH: CCGA Inner Results Bilinear Heuristic. 

This struct is made to store results from the bilinear reformulations using the FMPH (Binlinear search), 
which is used as another type of CCGAIR. 
"""
mutable struct IRBHeuristic <: CCGAIR
    "an instance of FMPH stepper that is used throughout the inner iterations of the CCGA. "
    fmphs::FMPHStepper 
    "The last FSP instance used. "
    fsp::FSP 
    "All the ds that had been passed to the FSP instance in the current iterations. "
    all_ds::Vector
    "All the q that had been passed to the FSP instance in the current iterations. "
    all_qs::Vector
    "All fmph objective values, forms one continuous trajectory. It's maybe not be monotonically decreasing. "
    upper_bounds::Vector
    """
    The objective of fsp, they form discrete trajectories. It's in the form of an array of array, inhomogenous. 
    Whenever FMPH failed to hit `ϵ` tolerance level for inner loop, a new empty vector will be inserted 
    to this lower bounds vector. 
    """
    lower_bounds::Vector{Vector}
    """
    A code representing different reasons that made the inner CCGA interminate: 
    * Status code  `1`: terminates because fsp is positive
    * status code  `0`: because fmp and fsp converged to zero. 
    * Status code `-1`: terminates because max iteration counter is reached. 
    """
    termination_status::Int 

    
    function IRBHeuristic()
        return new()
    end
end


"""
A function factored out and its being shared by instances of `IRBHeuristic`, and `IRBReform`. 
It reads out the values for the JuMP variables output from the FSP solve. 
"""
function ProduceReport(::CCGAIR, fsp::FSP)::Vector{String}
    # WARN:[?]() This code is not yet tested! 
    string_list = Vector{String}()
    u = fsp.u
    C = MatrixConstruct.C_
    # u = [c, p, h, g_plus, g_minus, dr, μ]
    var_coef_holder_list = MatrixConstruct.u
    for var_coef_holder in var_coef_holder_list
        (starting_at, ending_at) = MatrixConstruct.ColumnRegimeFor(C, var_coef_holder)
        for idx = starting_at:ending_at
            push!(string_list, "$(u[idx]) = $(u[idx]|>value)")
        end
    end
    return string_list
end

"""
    ProduceReport(this::IRBReform)::String

We report the feasible solutions produced by the FSP. Here is a list of parameters we are 
interested in: 
* `p, p',sr, sr',regu, regu', regd, regd' ,nsp, nsp',g+, g-, rho+, rho-`; 
And the values of these parameters will be returned as a formatted multi-line text. 
This function when the outer loop terminates and the inner loop had been identified 
as the last iteration of the inner loop. 

# Arguments
- `this::CCGAInnerResults`: It's a member method of this type. 
"""
function ProduceReport(this::IRBReform)::String
    # In addition to produce the report in text, these parameters need to be stored as flattend array after the algorithm is finished. 
    # This is for the CCGA outterloop.

    # WARN: [?]() This code is not yet tested. 
    string_list = Vector{String}()
    append!(string_list, ProduceReport(this, this.fsp))
    ρ⁺ = reshape(GetRhoPlus(this.fmp), size(MatrixConstruct.d))
    ρ⁻ = reshape(GetRhoMinus(this.fmp), size(MatrixConstruct.d))
    push!(string_list, "ρ⁺:")
    push!(string_list, repr("text/plain", ρ⁺))
    push!(string_list, "ρ⁻: ")
    push!(string_list, repr("text/plain", ρ⁻))
    # insert newline character text at the end of each line oof. 
    for idx in length(string_list) + 1 :-1: 1
        insert!(string_list, idx, "\n")
    end
    
return join(string_list) end

"""
    ProduceReport(this::IRBHeuristic)::String

We report the feasible solutions produced by the FSP. Here is a list of parameters we are 
interested in: 
* `p, p',sr, sr',regu, regu', regd, regd' ,nsp, nsp',g+, g-`; 
And the values of these parameters will be returned as a formatted multi-line text. 
The value `ρ⁺`, `ρ⁻` are not stored, nor is the last demand value tried by the instance of `FMPHStepper`. 

# Arguments
- `this::CCGAInnerResults`: It's a member method of this type. 
"""
function ProduceReport(this::IRBHeuristic)::String
    # TODO:[x](2) Produce Report for IRBHeuristic. 
    @warn "CCGAIRBH Producereport not yet implemented. "
    string_list = Vector{String}()
    append!(string_list, ProduceReport(this, this.fsp))
    # insert newline character text at the end of each line oof. 
    for idx in length(string_list) + 1 :-1: 1
        insert!(string_list, idx, "\n")
    end
    
return join(string_list) end




"""
Produce a figure that contains the objective values of the FMP and the FSP, plotted on the same graph. 
The x-axis is the number of iterations and the y-axis are the objective values of both FSP, FMP. 
The function will return a figure. 

# Arguments: 
- `this::CCGAInnerResults`: A member method of this type. 
"""
function ProducePlot(this::IRBReform)::Plots.Plot
    fig = plot(
        (1:length(this.upper_bounds)|>collect) .- 0.5, 
        this.upper_bounds;
        label="upper_fmp",
        marker=:x
    )
    plot!(
        fig, 
        1:length(this.lower_bounds)|>collect,
        this.lower_bounds;
        marker=:+,
        label="lower_fsp"
    )
return fig end

"""
Produce the plot for the inner iterations values for the fmphs. 
"""
function ProducePlot(this::IRBHeuristic)::Plots.Plot 
    #LATER: ProducePlot for the IRBHeuristic in the inner function.   
    @warn "Produce plot for CCGAIRBH not yet implemented. "
    return plot()
end

# ======================================================================================================================
# CCGA Outer Loop struct 
# ======================================================================================================================


"""
    CCGAOuterResults()

# CCGAOuterResults
Another struct that models the data exposed during the iterations of the full CCGA forloop (Inner and Outter). 
It will stores the following items: 
------
`inner_loops::Vector{CCGAInnerResults}`, `msp_objectives::Vector{Float64}`, `msp_gamma::Vector{Vector{Float64}}`, 
`msp::Union{MSP, Nothing}`, `termination_status::Int`, `fmp_initial_objectives::Vector{Float64}`, 
`fsp_initial_objectives::Vector{Floats}`
"""
mutable struct OutterResults
    
    """
    A list of CCGAInnerLoop that is obtained during the iterations of the outer CCGA loops. 
    They can be of different types depending on which function is used for the execution of the Inner CCGA. 
    """
    inner_loops::Vector{CCGAIR}
    "All the objective values of the msp for each iterations of the outter CCGA, initial msp objective without cut is in this vector. "
    msp_objectives::Vector{Float64}
    "The vector of gammas from the msp is stored here. "
    msp_gamma::Vector{Vector{Float64}}
    "Final instance of the master problem is stored here, it has all the cuts made for outer CCGA iterations. "
    msp::Union{MSP, Nothing}
    """
    A termination status is a code indicating the reasons of terminations for the outter loops of the algorithm. 
    - `-1`: max iteration is reached and it didn't break out of the forloop due to any conditions. 
    - `-2`: break out pre-maturally due to inner ccga forloop reaching max iterations counter. 
    - `0`: Successfully finished all iterations without reaching the outer maximum iterations. 
    """
    termination_status::Int
    "The initial objective value for the fmp at the start of each of the inner iterations. "
    fmp_initial_objectives::Vector{Float64}
    "The initial objective values for the fsp at the start of each of the inner CCGA for loops."
    fsp_initial_objectives::Vector{Float64}

    function OutterResults()
        this= new(
            Vector{IRBReform}(),
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
function (this::OutterResults)(that::IRBReform)::OutterResults
    push!(this.inner_loops, that)
    push!(this.fmp_initial_objectives, that.upper_bounds[1])
    push!(this.fsp_initial_objectives, that.lower_bounds[1])
return this end

"""
Pass an instance of the inner loop CCGA struct and transfer the results from the 
    inner iterations. More specifically the value of FSP and FMPH. 
"""
function (this::OutterResults)(that::IRBHeuristic)
    push!(this.inner_loops, that)
    push!(this.fmp_initial_objectives, that.upper_bounds[1])
    push!(this.fsp_initial_objectives, vcat(that.lower_bounds...)[end])
    return this 
end


"""
Call on an instance of `::MSP`, we will records all the gamma of the MSP instance, and then 
we will store the objective values of the current instance of **MSP**. 
"""
function (this::OutterResults)(that::MSP)::OutterResults
    push!(this.msp_gamma, GetGamma(that))
    push!(this.msp_objectives, objective_value(that))
return this end


"""

    ProduceReport(this::OutterResults)::String

Returns a *string* that is reporting all the results 
after the finishing the outer forloop iterations. We will be 
reporting the following important parameters: 
1. The initial objective values for the fmp, and fsp during each of the inner CCGA iterations. 
2. what are the upper bound and lower bound from the fmp and fsp during the last CCGA inner iteration. 
3. Make plots and return the figures if it's being asked to do it. 

### Arguments
- `this::CCGAOuterResults`: A member method for this type. 
"""
function ProduceReport(this::OutterResults)::String
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
    ProduceCSVFiles(this::OutterResults)

Print out all the decision variables for the last instance of the FSP, in a formatted CSV files. 
It will get the JuMP decision variable from an instance of `CCGAIR`. 

It will also print out all the objective values for the fsp at the last inner loop iterations as well. 
"""
function ProduceCSVFiles(this::OutterResults)::Nothing
    gammas_matrix = this.msp_gamma # this later. 
    u = this.inner_loops[end].fsp.u
    all_decision_var_vals = Vector{AbstractArray}()
    all_decision_var_names = Vector{Vector}()

    for var_coef_holder in MatrixConstruct.u
        (starting_at, ending_at) = MatrixConstruct.ColumnRegimeFor(MatrixConstruct.C_, var_coef_holder)
        decision_var_val = Vector{AbstractFloat}()
        for idx = starting_at: ending_at
            push!(decision_var_val, u[idx].|>value) # can you index the whole array using an array here? 
        end
        decision_var_val = reshape(decision_var_val, :, MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON)'
        decision_var_names = ["$(var_coef_holder.v)[$k]" for k in 1: size(decision_var_val, 2)]
        push!(all_decision_var_vals, decision_var_val)
        push!(all_decision_var_names, decision_var_names)
    end
    u_matrix = hcat(all_decision_var_vals...)
    u_matrix_col_attr = vcat(all_decision_var_names...)
    
    # Saving all the data into .csv files. 
    CSV.write("$(SESSION_DIR)/u_final.csv", Tables.table(u_matrix), header=u_matrix_col_attr)
    CSV.write("$(SESSION_DIR)/gammas_all.csv", Tables.table(hcat(gammas_matrix...)))
    CSV.write("$(SESSION_DIR)/gamma_objectives.csv",this.msp_objectives|>Tables.table)
return nothing end


"""
    SaveAllModels(this::CCGAOuterResults)

saves all the the models involved during the execution of the full CCGA 
iterations. The files will be all the final version of *FMP*, *MSP* for the inner and outer iterations in the `MOI`
file format compressed into `.gz`. 
"""
function SaveAllModels(this::OutterResults)
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
makes plots for `CCGAOutterResults`, there are 2 plots: 
1. All the final objective values of the FSP problems and the 
    final values of the FMP problems. 
2. The FSP and the MSP values for the last inner iterations of the CCGA.  
# Arguments: 
- `this::CCGAOuterResults`: It's a member method of this type. 
"""
function ProducePlots(this::OutterResults)
    fig1 = this.inner_loops[end]|>ProducePlot
    fig2 = plot(this.fmp_initial_objectives; label="fmp_initials_vals")
    plot!(fig2, this.fsp_initial_objectives; label="fsp_initial_vals")
    fig3 = plot(this.msp_objectives; label="msp_objectives")
return fig1, fig2, fig3 end


### ====================================================================================================================
### The actual Algorihm part for the CCGA. 
### ====================================================================================================================

"""
    CCGAInnerLoop(
        gamma_bar::Vector{N1}, 
        w_bar::Vector{N2},
        d_hat::Vector{N3};
        epsilon::Float64=0.1,
        max_itr::Int=8
    ) where {N1<:Number, N2<:Number, N3 <:Number}

This is the MIP reformulation of the Bi-linear problem. 

# Positional arguments: 
- `gamma_bar::Vector{N1}`: The initial gamma bound for each of the demand decision variable, predicted by the MSP.
- `w_bar::Vector{N2}`: The primary generator setup given by the *MSP*. 
- `d_hat::Vector{N3}`: The center of the demand interval.

# Named arguments: 
- `epsilon::Float64=0.1`: The absolute tolerance for terminating the inner CCGA for loop. 
- `max_iter::Int=8`: The maximum number of iterations before the for loop terminates. 
- `kwargs...`: This parameter is for ignoring invalid named parameters passed by other functions. 
"""
function InnerLoopMIP(
    gamma_bar::Vector{N1}, 
    w_bar::Vector{N2},
    d_hat::Vector{N3};
    inner_epsilon::Number=0.1,
    inner_max_itr::Int=8,
    kwargs...
) where {N1<:Number, N2<:Number, N3 <:Number}
    
    premise = "During executing CCGA Inner for loop: "
    @assert length(w_bar) == size(MatrixConstruct.B, 2) "$(premis)w̄ has the wrong size. please verify."
    @assert length(d_hat) == size(MatrixConstruct.H, 2) "$(premise)d̂ has the wrong size, please check the code. "
    @assert inner_epsilon >= 0 "$(premise)ϵ for terminating should be non-negative. "
    @assert inner_max_itr >= 0 "$(premise)Maximum iterations for the inner CCGA forloop should be non-negative. "
    @assert !(0 in (gamma_bar .>= 0)) "$(premise)γ̄ should be non-negative. "
    
    γ̄ = gamma_bar; d̂ = d_hat; ϵ=inner_epsilon; w̄ = w_bar
    # fill in these containers to return
    lowerbound_list = Vector{Float64}()
    upperbound_list = Vector{Float64}()
    all_qs = Vector{Vector}()
    all_ds = Vector{Vector}()
    termination_status = 0

    # ------------------------------------------------------------------------------------------------------------------
    # Preparing solvers
    model_fmp = GetJuMPModel(
        solver_name="FMP", 
        mip_focus=1, 
        log_to_console=1
    )
    AdaptSolverTimeout(model_fmp)
    fmp = FMP(w̄, γ̄, d̂, model_fmp)
    @info "$(TimeStamp()) Inner loop is initialized with fmp, and we are solving the initial fmp. "|>SESSION_FILE1
    Solve!(fmp)
    if objective_value(fmp)|>isnan
        DebugReport(fmp, "fmp_debug_report_inner_ccga")
        @assert false "$(premise) FMP is either unbounded or unfeasible on the first solve of FMP. "
    end
    push!(upperbound_list, objective_value(fmp))
    fsp = nothing
    for II in 1:inner_max_itr
        d = GetDemandVertex(fmp); push!(all_ds, d)
        model_fsp = GetJuMPModel()
        @info "$(TimeStamp()) FSP is made and we are solving it. "|>SESSION_FILE1
        fsp = FSP(w̄, d, model_fsp)
        Solve!(fsp); push!(lowerbound_list, fsp |> objective_value)
        @info "$(TimeStamp()) (FSP Lower, FMP Upper) = $((lowerbound_list[end], upperbound_list[end])) at itr = $II"|>SESSION_FILE1
        if lowerbound_list[end] > ϵ
            @info "$(TimeStamp()) Inner CCGA forloop terminated due to a positive lower of bound of"*
            ": $(lowerbound_list[end]) from FSP which is higher than: ϵ=$(ϵ)."|>SESSION_FILE1
            termination_status = 1
            break
        end
        q = Getq(fsp); push!(all_qs, q)
        IntroduceCut!(fmp, q)
        @info "$(TimeStamp()) CCGA Inner loop continues, constraints is introduced to fmp and we are solving it. "|>SESSION_FILE1
        AdaptSolverTimeout(model_fmp); Solve!(fmp)
        @assert !(objective_value(fmp)|>isnan) "$(premise) FMP"*
            " is infeasible or unbounded DURING the inner CCGA iterations. "
        push!(upperbound_list, objective_value(fmp))
        if abs(upperbound_list[end] - lowerbound_list[end]) < ϵ
            @info "$(TimeStamp()) Inner CCGA forloop termminated due to convergence of"*
                " FSP and FMP on tolerance level ϵ=$(ϵ), new fmp returns: $(fmp|>objective_value)"|>SESSION_FILE1
            break
        end
        
        if II == inner_max_itr
            termination_status = -1
        end
        
        model_fmp|>AdaptSolverTimeout

    end
    
return IRBReform(
    fmp,
    fsp,
    all_ds,
    all_qs,
    upperbound_list,
    lowerbound_list, 
    termination_status
) end


"""
This inner CCGA routine uses alternating heuristic to approximate the lower bound for the 
true value of `FMP`.

- `gamma_bar::Vector{N1}`: This is the uncertainty bounds returned by the master problem. 
- `w_bar::Vector{N2}`: Primary discrete variable from the master proble. 
- `d_hat::Vector{N3}`: The center of the uncertainty interval. 
- `epsilon::Float64=0.1`: The tolerance for deciding whether fmph is close to zero, or not close to zero. 
- `inner_max_itr::Int=8`: The maximal iterations allowed to compute the alternative heuristic for the FMPH. 
- `N::Int=10`: Controls the total number of attempts do a random search on demands for the FMPH system to bump up the 
heurstic each time when the value of heuristic estimate from fmph is lower than `epsilon`. 
- `M::Int=10`: It is the maximum number of times after a `q` is produced from FSP, and added as a cut for the FMPHs
and then performing the laternative heuristic based on the new `q`, as a results, it limits the maximum number of cuts 
can be added to the FMPHs, independent of the parameter `inner_max_itr`. 
- `kwargs...`: Leftover parameters for ignore extra named arguments passed from other functions. 

"""
function InnerLoopHeuristic(
    gamma_bar::Vector{N1},
    w_bar::Vector{N2},
    d_hat::Vector{N3};
    epsilon::Float64=0.1,
    inner_max_itr::Int=8, 
    N::Int=5,
    M::Int=10, 
    kwargs...
) where {N1<:Number, N2<:Number, N3 <:Number}
    premise = "During executing CCGA Inner for loop: "
    @assert length(w_bar) == size(MatrixConstruct.B, 2) "$(premis)w̄ has the wrong size. please verify."
    @assert length(d_hat) == size(MatrixConstruct.H, 2) "$(premise)d̂ has the wrong size, please check the code. "
    @assert epsilon >= 0 "$(premise)ϵ for terminating should be non-negative. "
    @assert inner_max_itr >= 0 "$(premise)Maximum iterations for the inner CCGA forloop should be non-negative. "
    @assert !(0 in (gamma_bar .>= 0)) "$(premise)γ̄ should be non-negative. "
    @info "$(TimeStamp()): Executing Inner loops with alternating bi-linear heuristic. "

    γ̄ = gamma_bar; d̂ = d_hat; ϵ = epsilon; w̄ = w_bar
    lowerbound_list = Vector{Vector{Float64}}() # for the fsp. 
    push!(lowerbound_list, Vector{Float64}())
    upperbound_list = Vector{Float64}() # for the fmph

    all_qs = Vector{Vector}() 
    all_ds = Vector{Vector}()
    local fmphs = FMPHStepper(w̄, γ̄, d̂, GetJuMPModel)
    push!(upperbound_list, fmphs|>objective_value)
    
    function AltUntilConverged(max_itr_alth::Int=20)
        previous = Inf; t = 0
        while (previous > fmphs|>objective_value) && (t < max_itr_alth)
            fmphs()
            previous = fmphs|>objective_value
            t += 1
        end
    end

    fsp = FSP(w̄::Vector{Float64}, fmphs|>GetDemands, GetJuMPModel())
    Solve!(fsp)  # always solved regardless of objective value of fmph. 
    push!(lowerbound_list, [fsp|>objective_value])
    m = 1
    termination_status = 0
    break_flag = false
    for II in 1:inner_max_itr
        if fmphs|>objective_value > ϵ
            if fsp|>objective_value > ϵ
                @info "$(TimeStamp()): FSP objective value: $(fsp|>objective_value) > $ϵ;\n"*
                "FMPHS obj value: $(fmphs|>objective_value) > $ϵ Terminate. "
                termination_status = 1  # exits due to FSP > ϵ
                break_flag = true
            else
                m += 1
                fsp = FSP(w̄::Vector{Float64}, fmphs|>GetDemands, GetJuMPModel()); Solve!(fsp)
                fmphs(fsp|>Getq)
                @info "$(TimeStamp()): FSP objective value $(fsp|>objective_value) < $ϵ, FMPHS obj > $(ϵ). Introducing cut to fmph and perform alt heuristic. "
                AltUntilConverged()
                if m >= M
                    @info "$(TimeStamp()): Terminates due to m=$m>=$M. "
                    termination_status = -1  # stagnations. 
                    break_flag = true
                end
            end
            # updates the lower bound and the upper bounds. 
            push!(upperbound_list, fmphs|>objective_value)
        else
            @info "$(TimeStamp()): FMPH has objective: $(fmphs|>objective_value), too low, we will try new demands to bump it up. "
            for _ = 1:N
                AltUntilConverged()
                if fmphs|>objective_value < ϵ # The suspicious case. 
                    @info "$(TimeStamp()): FMPH objective: $(fmphs|>objective_value), trying new random demands. "
                    fmphs|>TryNewDemand
                else
                    @info "FMPH objective: $(fmphs|>objective_value) exceed ϵ, done and exit inner heuristic loop."
                    break # break out of this forloop and continue in the outer forloop. 
                   
                end
            end
            push!(lowerbound_list, Vector{Float64}())  # add an empty vector to the lower bounds produced by the FSP instance. 
        end
        push!(upperbound_list, fmphs|>objective_value)
        # before we break... 
        push!(all_ds, fmphs|>GetDemands)
        push!(all_qs, fsp|>Getq)
        if break_flag
            break
        end
    end
    # if forloop exists with max_itr == II, it counts as success. 
    results = IRBHeuristic()
    results.termination_status = termination_status
    results.fsp = fsp
    results.fmphs = fmphs
    results.lower_bounds = lowerbound_list
    results.upper_bounds = upperbound_list
    results.all_ds = all_ds
    results.all_qs = all_qs

return results end

"""
    RoutineFor(msp::MSP, inner_loop_results::IRBHeuristic)::MSP

A subroutine that introduces the cut to the master problem given that the subroutine for the inner 
loop computes using the alternating search heuristic. 
"""
function RoutineFor(msp::MSP, inner_loop_results::IRBHeuristic)::MSP
    return IntroduceCut!(msp, inner_loop_results.fmphs|>GetDemands)
end

"""
    RoutineFor(msp::MSP, inner_loop_results::IRBReform)::MSP

A subroutine that itroduces the cut to the master problem given that the subroutine for the inner loop 
computes using the bilinear reformulation with MIP
"""
function RoutineFor(msp::MSP, inner_loop_results::IRBReform)::MSP
    IntroduceCut!(
        msp, 
        GetRhoPlus(inner_loop_results.fmp), 
        GetRhoMinus(inner_loop_results.fmp), 
    )
    return msp
end



"""

    function CCGAOuterLoop(
        d_hat::Vector{N1}, 
        gamma_upper::N2;
        outter_epsilon::N4=0.1,
        outer_max_itr::Int=40,
        make_plot::Bool=true, 
        inner_routine::Function=CCGAInnerLoop,
        kwargs...
    ) where {N1 <: Number, N2 <: Number, N4<:Number}

Performs the outter forloop of the CCGA algorithm with initialized parameters.

# Positional Arguments
- `d_hat::Vector{N1}`: The center of the demands uncertainty interval. 
- `gamma_upper::N2`: The initial scalar upper bound for all the γ in the uncertainty interval. 

# Keyword Arguments:
- `inner_epsilon`: The tolerance for the lower bound and upper bound between FMP, FSP, and it's used for the termination 
    conditions for the CCGA inner forloop. 
- `make_plot::Bool=true`: Whether to make a plot for all the results obtain from the execution of the inner CCGA forloop. 
- `inner_max_itr=8`: The maximum time s for executing the inner forloop of CCGA. 
- `outer_max_itr=40`: The maximum number of tiles for excuting the out forloop of CCGA. 
- `make_plot::Bool=true`: Make plots for each of the inner CCGA iterations while running, also make plots after 
    everything is finished. 
- `objective_types::Int=0`: See doc for *MSP* for more. 
- `block_demands::Int=1`: See doc for type *MSP* for more information. 
- `inner_routine::Function=CCGAInnerLoop`, the inner routine function you want the outter CCGA forloop function 
to call. 
- `session_time_out::Bool=False`
"""
function OuterLoop(
    d_hat::Vector{N1}, 
    gamma_upper::N2;
    outer_max_itr::Int=40,
    make_plot::Bool=true, 
    inner_routine::Function=InnerLoopMIP,
    session_time_out::Int=-1,
    kwargs...
) where {N1 <: Number, N2 <: Number}
    context = "During the execution of the outter loop of CCGA: "
    @assert length(d_hat) == size(MatrixConstruct.H, 2) "$(context)"*
    "The length of d_hat is wrong, it's $(length(d_hat)), but we want: $(size(MatrixConstruct.H, 2)). "
    @assert gamma_upper > 0 "$context gamma_upper should be a strictly positive number."
    # Store the parameters of the CCGA for reports. 
    let
        # kwargsd = Dict([(item[1], item[2]) for item in kwargs])
        SESSION_FILE2() do io
            # Problem parameters for algorithm 
            println(io, 
                foldl(*, [("$(item[1])=$(item[2])\n") for item in kwargs])
            )
            # Problem parameters for the systems
            println(io, "Gamma upper, or M is: $(gamma_upper)")
            println(io, "d̂: ")
            println(io, repr("text/plain", d_hat))
            println(io, "HORIZON = $(MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON)")
            println(io, "budget = $(MatrixConstruct.CONST_PROBLEM_PARAMETERS.Φ)")
            println(io, "H̄ = $(MatrixConstruct.STORAGE_SYSTEM.Capacity)")
            println(io, "Pmax = $(MatrixConstruct.PRIMARY_GENERATORS.Pmax)")
            println(io, "Pmin = $(MatrixConstruct.PRIMARY_GENERATORS.Pmin)")
            println(io, "RU = $(MatrixConstruct.PRIMARY_GENERATORS.RU)")
            println(io, "RD = $(MatrixConstruct.PRIMARY_GENERATORS.RD)")
            println(io, "DR = $(MatrixConstruct.DEMAND_RESPONSE.R)")
        end
    end
    if session_time_out >= 0
        global SESSION_START_TIME = datetime2unix(now())|>floor|>Int
        global SESSION_TIME_OUT = session_time_out
        @info "Session time limit has been set to $(SESSION_TIME_OUT) seconds. "*
            "This will dynamically adjust solver time_out options depending on the amount of time remainds. "
    end
    γ⁺ = gamma_upper; d̂ = d_hat
    model_mp = GetJuMPModel()
    mp = MP(model_mp, γ⁺)
    model_msp = GetJuMPModel(solver_name="MSP", log_to_console=1, optimality_gap=0.001)
    AdaptSolverTimeout(model_msp)
    msp = MSP(model_msp, d̂, γ⁺; kwargs...)
    PortOutVariable!(mp, :d) do d fix.(d, d̂, force=true) end
    PortOutVariable!(mp, :v) do v fix.(v, 0, force=true) end
    Solve!(mp)
    w̄ = Getw(mp)
    Solve!(msp)
    γ̄ = GetGamma(msp)
    outer_counter = 0
    outer_results = OutterResults()
    for _ in 1:outer_max_itr
        outer_counter += 1
        @info "$(TimeStamp()) Outter Forloop itr=$outer_counter" |> SESSION_FILE1
        Results = inner_routine(γ̄, w̄, d̂; kwargs...)
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
        if Results.termination_status == 0
            @info "$(TimeStamp()) Outer for loop terminated due to convergence of FMP, FSP to an objective value of zero."|>SESSION_FILE1
            break
        end
        @info "$(TimeStamp()) Introduced cut to the msp and we are solving it. "|>SESSION_FILE1
        RoutineFor(msp, Results) # Deploy the cut routine based on the inner results type, for the master problem. 
        AdaptSolverTimeout(model_msp)
        Solve!(msp)
        outer_results(msp)
        w̄ = Getw(msp)
        γ̄ = GetGamma(msp)
        @info "$(TimeStamp()) Objective value of msp settled at: $(objective_value(msp)). "|>SESSION_FILE1
        @assert !isnan(objective_value(msp)) "$context objective value for msp is NaN. "
        @assert !isinf(objective_value(msp)) "$contex objective value of msp is inf. "
        outer_results|>ProduceCSVFiles
    end

    outer_results.msp = msp
    if outer_counter == outer_max_itr 
        outer_results.termination_status = -1
    end

    # Print the results to files! 
    SESSION_FILE3() do io
        write(io, outer_results|>ProduceReport)  # print out master problem report first! 
        write(io, outer_results.inner_loops[end]|>ProduceReport) # print out the inner loop results next. 
    end
    

    if inner_routine === InnerLoopMIP
        SaveAllModels(outer_results)
    end

    if make_plot
        figs = outer_results|>ProducePlots  
        savefig(figs[1], SESSION_DIR*"/"*"last_fmp_fsp")
        savefig(figs[2], SESSION_DIR*"/"*"initial_fmp_fsp")
        savefig(figs[3], SESSION_DIR*"/"*"objectives_msp")
    end
    
    # END

return outer_results end




