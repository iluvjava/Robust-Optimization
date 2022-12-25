### ====================================================================================================================
### Main Problem:
###   * Search for primary feasible solution given branch cut and bounds.
###   * Is able to perform a feasibility search for all the involved variables in the model. Both primary and secondary
###     variables.
###     * Is able to produce report.
### ====================================================================================================================

"""
MP: Master problem. 

# Constructor 
    MP(M::Model, gamma_upper=1e4)

"""
@ProblemTemplate mutable struct MP <: Problem
    
    w::Array{VariableRef, 3}
    gamma::Vector{VariableRef}
    "an upper bound for the gamma variable: γ⁺"
    gamma_upper::Number         
    "secondary continuous decision variables: u"
    u::Vector{VariableRef}      
    "secondary discrete decision variables."
    q::Vector{VariableRef}
    "the demand variables, continuous in the case of master problem."
    d::Vector{VariableRef}      
    "The slack decision variable."
    v::Vector{VariableRef} 
    "Given number of generators. "
    G::Int64
    "The time horizon parameter. "
    T::Int64
    "Min down time for primary generators. "
    Tmind::Array{Int}           
    "Min up time for the primary generators. "
    Tminu::Array{Int}

    # settings and options stuff:
    Verbose::Bool

    function MP(M::Model, gamma_upper=1e4)
        this = new(M)
        this.G = MatrixConstruct.PRIMARY_GENERATORS.generator_count
        this.T = MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON
        this.Tmind = MatrixConstruct.PRIMARY_GENERATORS.Tmind
        @assert any(this.T .> this.Tmind) "Tmind: Minimum down time has to be less than "*
        "time horizon"
        this.Tminu = MatrixConstruct.PRIMARY_GENERATORS.Tminu
        @assert any(this.T .> this.Tminu) "Tmind: Minimum up time has to be less than "*
        "time horizon"
        this.gamma_upper = gamma_upper
        this.con = Vector()
        EstablishMainProblem!(this)
    return this end

    function MP(gamma_upper=1e4) return MP(Model(HiGHS.Optimizer), gamma_upper) end

end


"""
    IntroduceVariables!(this::MP)

Introduce variables to the MP's model:
* Secondary continuous decision variables: u
* Secondary discrete decision variables: q
* demands variables: d
"""
function IntroduceVariables!(this::MP)
    model = this |> GetModel
    this.u = PrepareVariablesForTheModel!(model, :u)
    this.q = PrepareVariablesForTheModel!(model, :q)
    this.d = @variable(model, d[1:size(MatrixConstruct.H, 2)] >= 0)
    this.v = @variable(model, v[1:size(MatrixConstruct.H, 1)] >= 0)
return this end



"""
    AddMainProblemConstraints!(this::MP)

Craete a constraints to check whether there exists an initial feasiblity solutions to the system.
"""
function AddMainProblemConstraints!(this::MP)
    model = this|>GetModel
    w = Flattenw(this)
    u = this.u
    q = this.q
    d = this.d
    B = MatrixConstruct.B
    C = MatrixConstruct.C
    G = MatrixConstruct.G
    H = MatrixConstruct.H
    h = MatrixConstruct.h
    γ = this.gamma  # Unused. 
    v = this.v
    push!(this.con, @constraint(model, B*w + C*u + G*q + H*d - v .<= h)...)
    # push!(this.con, @constraint(model, d .>= γ)...)
return this end


"""
    DemandFeasible(this::MP, demand::Vector{N}) where {N<:Number}

Test if a given demand vector is feasible for the main problem.
* Set the demands
* Solve
* Get results
* Unfix the demends vector.
"""
function DemandFeasible(this::MP, demand::Vector{N}) where {N<:Number}
    d = this.d
    fix.(d, demand; force=true)
    @objective(this.model, Min, sum(this.v))
    Solve!(this)
    ToReturn = !((this |>objective_value) > 0 || (this|>objective_value|>isnan))
    unfix.(d)
return ToReturn end

"""
    DemandFeasible(this::MP, demand::N)  where {N<:Number}

Adapting the function when the user decided to put into a number instead 
of a vector for the demand to test. 
"""
function DemandFeasible(this::MP, demand::N)  where {N<:Number}
return DemandFeasible(this, demand*ones(size(this.d))) end



"""
Creates the main problem objectives:
    * Find a maximum lower bound for demands on all buses, all time.
    * Fix the slack variable to be zero.
"""
function MainProblemObjective!(this::MP)
    @warn "METHOD DEPRECATED. "
    m = this|>GetModel
    γ = m[:γ]
    v = this.v
    d = this.d
    # @objective(m, Max, γ)
    # fix.(v, 0; force=true)
return end


"""
Establish the main problem, and then solving it will provide some suggestions for the
best average demands for the system.
* The constraints
* The objective
* The variables as well.
"""
function EstablishMainProblem!(this::MP)
    # Prepre Variables for the main problem.
    IntroduceMasterProblemVariables!(this)
    IntroduceVariables!(this)

    # Prepare constraints for the Main Problem.
    PreppareConstraintsPrimary!(this)
    AddMainProblemConstraints!(this)

    # Add it yourself whenever you make it. 
    # MainProblemObjective!(this)
return this end



### ====================================================================================================================
#   Master Problem:
#       * Determine primary decision variables only.
#       * Similar format compare to the Main Problem.
# Options and settings
#   objective_types: 
#       1. maximize gamma min.
#       2. maximize sum of all gamma.
#   block_demands_types
#       0. no restriction on demand interval variable gamma. 
#       1. restrict all generator at different time to have the same demand intervals. 
### ====================================================================================================================

@ProblemTemplate mutable struct MSP <: Problem

    w::Array{VariableRef, 3}
    d_hat::Array{Float64}
    gamma::Vector{VariableRef}
    gamma_upper::Number         # an upper bound for the gamma variable.
    gamma_min::VariableRef       # A minimum bound for gamma. 
    G::Int64; T::Int64          # Given number of generators and time horizon for the problem.
    Tmind::Array{Int}           # Min up down for primary generators
    Tminu::Array{Int}           # Min up time for primary generators

    cut_count::Int
    u::Vector{Vector{VariableRef}}  # Decision variables from the cut introduced to the master problem. 
    q::Vector{Vector{VariableRef}}
    # artificial_bound_at::Int        # an artificla bounds for the objective, it's undefined before adding any cuts 
    
    block_demands_types::Int    # whether to lock the group of demand uncertainty interval gamma to the same value. 
    objective_types::Int            # a integer that switches between different type of objectives for the msp problem. 

    function MSP(
        model::Model, 
        d_hat::Vector{T}, 
        gamma_upper;
        block_demands::Int = 0, 
        objective_types::Int = 1
    ) where {T<:AbstractFloat}
        this = new()
        this.model= model
        this.gamma_upper = gamma_upper
        this.G = MatrixConstruct.PRIMARY_GENERATORS.generator_count
        this.T = MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON
        this.Tmind = MatrixConstruct.PRIMARY_GENERATORS.Tmind
        @assert any(this.T .> this.Tmind) "Tmind: Minimum down time has to be less than "*
        "time horizon"
        this.Tminu = MatrixConstruct.PRIMARY_GENERATORS.Tminu
        @assert any(this.T .> this.Tminu) "Tmind: Minimum up time has to be less than "*
        "time horizon"
        this.gamma_upper = gamma_upper
        this.con = Vector()
        this.d_hat = d_hat
        this.cut_count = 0
        this.u = Vector{Vector{VariableRef}}()
        this.q = Vector{Vector{VariableRef}}()
        this.block_demands_types = block_demands;
        this.objective_types = objective_types;
    
        IntroduceMasterProblemVariables!(this)
        PreppareConstraintsPrimary!(this)
        MasterProblemObjective!(this)

    return this end

    function MSP(d_hat::Vector{T}, gamma_upper=50) where {T<:AbstractFloat}
    return MSP(Model(HiGHS.Optimizer), d_hat, gamma_upper) end

end


"""
Create the upper bound for demands for the master problem.
    * Adds the demands bounds constraints
    * Adds the objective bounds maximizations objectives for the model.
"""
function MasterProblemObjective!(this::MSP)
    model = GetModel(this)
    γ = this.gamma; γmin = this.gamma_min

    if this.objective_types == 1
        @objective(model, Max, γmin)
    elseif this.objective_types == 2
        @objective(model, Max, sum(γ))
    else
        error("this.objective_types = $(this.objective_types) is not a valid argument. ")
    end
return this end


# METHOD SHARING WITH OTHER TYPES --------------------------------------------------------------------------------------


"""
    IntroduceMasterProblemVariables!(this::Union{MP, MSP})

Introduce the binary decision variable w for the primary generator. For instance of both `MP`, `MSP`, it constructs: 
* w
* γ
* γmin, a lower bound for block/max min type of objective for the MSP problem. 
"""
function IntroduceMasterProblemVariables!(this::Union{MP, MSP})
    model = this |> GetModel
    this.w = @variable(
        model,
        w[
            [:x, :y, :z],
            1:this.G,
            1:this.T
        ], Bin
    )

    # Scalar Continuous variables for demand interval.
    time_horizon = MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON
    # B̄ = size(MatrixConstruct.H, 2)
    this.gamma = @variable(
        model,
        γ[1:MatrixConstruct.B̄, 1:time_horizon],
        lower_bound=0, upper_bound=this.gamma_upper
    )[:]
    
    if isa(this, MSP)
        this.gamma_min = @variable(
            model, 
            γmin,
            lower_bound = 0
        )
    end

return this end



"""
Prepare the constraints for the primary generator, the system
    * Aw <= b
"""
function PreppareConstraintsPrimary!(this::Union{MP, MSP})
    model = this|>GetModel
    w = (model)[:w]

    push!(
        this.con,
        @constraint(
            this|>GetModel, w[:x, :, :] + w[:z, :, :] .<= 1
        )...
    )

    push!(
        this.con,
        @constraint(this|>GetModel,
            w[:y, :, 1] .== w[:x, :, 1] - w[:z, :, 1]
        )...
    )

    for t in 2: this.T
        push!(
            this.con,
            @constraint(
                this|>GetModel,
                w[:y, :, t] - w[:y, :, t - 1] .== w[:x, :, t] - w[:z, :, t]
            )...
        )
    end

    for n in 1:this.G, t in this.Tminu[n]: this.T
        push!(
            this.con,
            @constraint(
                this|>GetModel,
                sum(
                    w[:x, n, tau] for tau in t - this.Tminu[n] + 1: t
                )
                <=
                w[:y, n, t]
            )
        )
    end

    for n in 1:this.G, t in this.Tmind[n]: this.T
        push!(
            this.con,
            @constraint(
                this|>GetModel,
                sum(
                    w[:z, n, tau] for tau in t - this.Tmind[n] + 1: t
                )
                <=
                1 - w[:y, n, t]
            )
        )
    end

    if isa(this, MSP)
        # depends on the option settings, we assert different type of block constraints on the demand invervals decision 
        # variable gamma. 
        if this.block_demands_types == 0
            # do nothing here. 
        elseif this.block_demands_types == 1
            num_gen = MatrixConstruct.B̄; t_horizon = MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON
            γ = this.gamma
            for t in 0:t_horizon - 1
                γ_block = γ[num_gen*t + 1: num_gen*(t + 1)]
                for (γ1, γ2) in zip(γ_block[1:end - 1], γ_block[2: end])
                    push!(this.con, @constraint(model, γ1 == γ2))
                end
            end
        else
            error("this.block_demands_types=$(this.block_demands_types) is an invalid argument. ")
        end

        # this constraint has to be the LAST constraints to be added to avoid any type of errors. 
        # this.artificial_bound_at = length(this.con) + 1   # an index to locate the artificial bound introduced. 
        # push!(this.con, @constraint(model, this.gamma_min .<= this.gamma)...)
    end
return this end



"""
    IntroduceCut!(
        this::MSP,
        rho_plus::Vector{Float64},
        rho_minus::Vector{Float64},
    )

Introduce feasibility cut for the master problem from the CCGA results. 

# Arguments
- `this::MSP`
- `rho_plus::Vector{float64}`: The rho plus solved from the FMP. 
- `rho_minus::Vector{Float64}`: The rho minus from the FMP. 

"""
function IntroduceCut!(
    this::MSP,
    rho_plus::Vector{Float64},
    rho_minus::Vector{Float64},
)
    model = this|>GetModel
    w = this |> Flattenw
    B = MatrixConstruct.B
    h = MatrixConstruct.h
    G = MatrixConstruct.G
    C = MatrixConstruct.C
    H = MatrixConstruct.H
    d̂ = this.d_hat
    ρ⁺ = rho_plus
    ρ⁻ = rho_minus
    γ = this.gamma

    Γ = γ|>Diagonal             
    this.cut_count += 1
    
    u = PrepareVariablesForTheModel!(model, :u, this.cut_count)
    q = PrepareVariablesForTheModel!(model, :q, this.cut_count)
    push!(this.u, u)
    push!(this.q, q)

    CutConstraints = @constraint(
        model, 
        B*w + C*u + G*q + H*d̂ + H*Γ*(ρ⁺ - ρ⁻).<= h, 
        base_name="Cut $(this.cut_count)"
    )

    push!(
        this.con,
        CutConstraints...   
    )
    
return CutConstraints end


# """
# Delete all the previous cut introduced to the master problem. 
# * There will be an error if no cut is introduced and we tried to delete them. 
# """
# function DeleteAllPreviousCut!(this::MSP)
#     startingAt = this.artificial_bound_at + 1
#     endingAt = length(this.con)
#     delete.(this.model, this.con[startingAt: endingAt])
#     deleteat!(this.con, startingAt: endingAt)
# return this end


"""
    function Getw(this::Union{MP, MSP})

Get the primary discrete decision variable as a vector of numbers for
the CCGA algorithm.
* It returns vector as decision variable
* The vector is w flattend into x, y, z, column major flattening.
* The returned type is Vector{VariableRef}
"""
function Getw(this::Union{MP, MSP})
return this|>Flattenw.|>value end

"""
    GetGamma(this::Union{MP, MSP})

Get the gamma vector from the master problem that is here. The gamma vector will be a float64 vector 
with the same length as the number of time horizon that is there. 
"""
function GetGamma(this::Union{MP, MSP})
return this.gamma.|>value end

"""
    Flattenw(this::Union{MP, MSP})

Fletten the w decision variable. 
"""
function Flattenw(this::Union{MP, MSP})
    Vec = Vector{JuMP.VariableRef}()
    push!(Vec, this.w[1, :, :][:]...)
    push!(Vec, this.w[2, :, :][:]...)
    push!(Vec, this.w[3, :, :][:]...)
return Vec end