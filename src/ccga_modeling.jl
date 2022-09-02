### ====================================================================================================================
### Important functions for variables constructions for JuMP models.
### ====================================================================================================================

using Mixers, Kronecker

"""
    Does a cartesian outter product on a list of vectors
    and then return all the products in the form of a tuple packed into
    a long long list.
"""
function CartesianOutterProductList(v_list::AbstractVector...)
    result = Iterators.product(v_list...)|>collect
return result[:] end


"""
    Given a variable coefficient holder, using its dimension to convert
    it to a list of tuples representing all the possible indexing of this
    multi-dimensional variable.
        * k: The extra forloop indexing from the CCGA algorithm.

    Returns:
        All possible indices for the decision variables represented by the
        instance of the coefficient holder.
"""
function IndicesList(
    holder::MatrixConstruct.VariableCoefficientHolder,
    k::Union{Nothing, Int}=nothing
)
    NdimList = [1:III for III in size(holder)]
    if k === nothing
        return CartesianOutterProductList(NdimList...)
    end
    push!(NdimList, k:k)
return CartesianOutterProductList(NdimList...) end


abstract type Problem
    # has a JuMP model in it.
    # MUST IMPLEMENT: GetModel(::Problem)
end

@premix mutable struct ProblemTemplate
    model::Model
    con::Vector{ConstraintRef}
end



"""
    Given a JuMP model, prepare the q, u variable for that model.
        * Returns the variable 'u' or 'q' packed into Vector{JuMP.VariableRef}.
"""
function PrepareVariablesForTheModel!(
    model::Model,
    variable::Symbol,
    ccga_itr::Union{Nothing, Int}=0
)
    ccga_itr = ccga_itr == 0 ? "" : "[$(ccga_itr)]"
    if variable == :u
        u = Vector{JuMP.VariableRef}()
        for v in MatrixConstruct.u
            push!(u, @variable(model, [IndicesList(v)], base_name="$(v.v)$(ccga_itr)", lower_bound=0)...)

        end
        return u
    end
    
    if variable == :q
        q = Vector{JuMP.VariableRef}()
        for v in MatrixConstruct.q
            push!(q, @variable(model, [IndicesList(v)], Bin, base_name="$(v.v)$(ccga_itr)")...)

        end
        return q
    end

return nothing end



### ====================================================================================================================
### Main Problem:
###   * Search for primary feasible solution given branch cut and bounds.
###   * Is able to perform a feasibility search for all the involved variables in the model. Both primary and secondary
###     variables.
###     * Is able to produce report.
### ====================================================================================================================

@ProblemTemplate mutable struct MP <: Problem
    
    w::Array{VariableRef, 3}
    gamma::Vector{VariableRef}
    gamma_upper::Number         # an upper bound for the gamma variable.
    u::Vector{VariableRef}      # secondary continuous decision variables.
    q::Vector{VariableRef}      # secondary discrete decision variables.
    d::Vector{VariableRef}      # the demand variables, continuous in the case of master problem.
    v::Vector{VariableRef}      # The slack decision variable.
    
    G::Int64; T::Int64          # Given number of generators and time horizon for the problem.
    Tmind::Array{Int}           # Min up down for primary generators
    Tminu::Array{Int}           # Min up time for primary generators

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
    Introduce variables to the model:
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
### ====================================================================================================================

@ProblemTemplate mutable struct MSP <: Problem

    w::Array{VariableRef, 3}
    d_hat::Array{Float64}
    gamma::Vector{VariableRef}
    gamma_upper::Number         # an upper bound for the gamma variable.
    G::Int64; T::Int64          # Given number of generators and time horizon for the problem.
    Tmind::Array{Int}           # Min up down for primary generators
    Tminu::Array{Int}           # Min up time for primary generators

    cut_count::Int
    u::Vector{Vector{VariableRef}}  # Decision variables from the cut introduced to the master problem. 
    q::Vector{Vector{VariableRef}}
    artificial_bound_at::Int         

    function MSP(
        model::Model, 
        d_hat::Vector{T}, 
        gamma_upper;
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
    γ = this.gamma
    @objective(model, Max, sum(γ))
return this end


# METHOD SHARING WITH OTHER TYPES --------------------------------------------------------------------------------------


"""
    Introduce the binary decision variable w for the primary generator.
        * w
        * γ
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
    #TODO: Change here for gamma
    this.gamma = @variable(
        model,
        γ[1:time_horizon],
        lower_bound=0, upper_bound=this.gamma_upper
    )[:]

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
        this.artificial_bound_at = length(this.con) + 1
    end
return this end



"""
    Introduce feasibility cut for the master problem which comes from the CCGA results.
        * Delete all constraints established for the MainProblem.
        * Continuous decision variable from FSP: u
        * Discrete decision variable from FSP: q
        * Adversarial demands from FMP: ρ⁺, ρ⁻
        * introduce new objective for maximum demands intervals satisfactions.
"""
function IntroduceCut!(
    this::MSP,
    rho_plus::Vector{Float64},
    rho_minus::Vector{Float64},
    aritifical_bound::Float64
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
    Σγ⁺ = aritifical_bound

    Γ = kronecker(γ|>Diagonal, ones(MatrixConstruct.B̄)|>Diagonal)|>collect|>Diagonal #TODO: Change here for gamma
    this.cut_count += 1
    
    u = PrepareVariablesForTheModel!(model, :u, this.cut_count)
    q = PrepareVariablesForTheModel!(model, :q, this.cut_count)
    push!(this.u, u)
    push!(this.q, q)

    # if v === nothing
    #     v = zeros(length(h))
    # end
    # model[:s] = s = @variable(
    #     model, 
    #     [1:length(h)], 
    #     lower_bound=0, 
    #     base_name="s"
    # )
    # fix.(model[:s][374:end], 0, force=true)
    # fix.(model[:s][1:381], 0, force=true)
    # Infiltrator.@exfiltrate
    
    if this.cut_count > 1
        previuosArtificialBound = popat!(this.con, this.artificial_bound_at)
        delete(model, previuosArtificialBound)
    end

    insert!(this.con, 
        this.artificial_bound_at, @constraint(model, sum(γ) <= Σγ⁺)
    )

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

"""
    Delete all the previous cut introduced to the master problem. 
    * There will be an error if no cut is introduced and we tried to delete them. 
"""
function DeleteAllPreviousCut!(this::MSP)
    startingAt = this.artificial_bound_at + 1
    endingAt = length(this.con)
    delete.(this.model, this.con[startingAt: endingAt])
    deleteat!(this.con, startingAt: endingAt)
return this end


"""
    Get the primary discrete decision variable as a vector of numbers for
        the CCGA algorithm.
    * It returns vector as decision variable
    * The vector is w flattend into x, y, z, column major flattening.
    * The returned type is Vector{VariableRef}

"""
function Getw(this::Union{MP, MSP})
return this|>Flattenw.|>value end

"""
    Get the gamma vector from the master problem that is here. The gamma vector will be a float64 vector 
    with the same length as the number of time horizon that is there. 
"""
function GetGamma(this::Union{MP, MSP})
return this.gamma.|>value end

"""
    Fletten the w decision variable. 
"""
function Flattenw(this::Union{MP, MSP})
    Vec = Vector{JuMP.VariableRef}()
    push!(Vec, this.w[1, :, :][:]...)
    push!(Vec, this.w[2, :, :][:]...)
    push!(Vec, this.w[3, :, :][:]...)
return Vec end




### ====================================================================================================================
# FSP: Lower bound searcher!
#   * Takes in a demands and it tries to satisfies it by choosing the secondary.
#   * discrete and contiuous variables.
### ====================================================================================================================
"""
    Given demand from FMP, test how feasible it is to determine an Lower Bound
    for the feasibility slack variables.
"""
@ProblemTemplate mutable struct FSP <: Problem
    # model::Model
    u::Vector{VariableRef}
    v::Vector{VariableRef}
    q::Vector{VariableRef}

    # Given parameters
    w::Vector{Float64}
    d_star::Vector{Float64} 

    sparse_vee::Bool

    function FSP()
        @warn("This construction only exists for testing purposes!")
        return FSP(
            zeros(size(MatrixConstruct.B, 2)),
            10,
            ones(MatrixConstruct.d|>length)
        )
    end

    """
        Constrct the FSP. 
        Parameter: 
            * w̄::Vector{Float64}
            * γ̄::Number
            * d⋆::Vector{Float64}
            * model::Model=Model(HiGHS.Optimizer)
        Parameters called by names: 
            * sparse_vee::Bool=false
    """
    function FSP(
        w::Vector{Float64}, 
       #  gamma::Vector{Float64},
        d_star::Vector{Float64}, 
        model::Model=Model(HiGHS.Optimizer);
        sparse_vee::Bool=false
    )
        this = new()
        this.model = model
        # this.gamma = gamma
        this.w = w
        this.d_star = d_star
        this.con = Vector{JuMP.ConstraintRef}()
        this.sparse_vee = sparse_vee
        this|>IntroduceVariables!
        this|>AddConstraints!
        @objective(this.model, Min, sum(this.v))
    return this end

    """
        Prepare all the variables for the FSP model, decision variables:
            * u
            * v
            * q
        As a list of flattened JuMP variable refs, it will store the flatten
        variables to the field of this struct.
    """
    function IntroduceVariables!(this::FSP)
        this.u = PrepareVariablesForTheModel!(this|>GetModel, :u)
        this.q = PrepareVariablesForTheModel!(this|>GetModel, :q)
        this.v = @variable(this.model, v[1:length(MatrixConstruct.h)], lower_bound=0)
        starting, ending = MatrixConstruct.CON_GROUPS["Demand Balance"]
        if this.sparse_vee
            set_upper_bound.(
                v[
                    setdiff(
                        1:size(MatrixConstruct.H, 1)|>Set,
                        starting:ending|>Set
                    )|>collect|>sort
                ], 
                0
            )
        end
    return end


    function AddConstraints!(this::FSP)
        u = this.u
        v = this.v
        q = this.q
        h = MatrixConstruct.h
        d = this.d_star
        w = this.w
        C = MatrixConstruct.C
        H = MatrixConstruct.H
        B = MatrixConstruct.B
        G = MatrixConstruct.G
        # @info "Preparing constraints for the FSP model. "
        push!(
            this.con, 
            @constraint(this.model, C*u + H*d + B*w + G*q - v .<= h)...
        )
    return end


    function Base.getindex(this::FSP, index::Symbol)
    return this.model[index] end

end


# ======================================================================================================================
# FMP, the Upper bound locator.
#   * Obtains a upper bound by giving a demands that can break the feasibility for all discrete secondary
#     decision variables.
# Supports:
#   * Solving the system to obtain
# ======================================================================================================================

abstract type AbsFMP <: Problem
    # FMP like problem. Different way of working with solving the FMP. 
end

@premix mutable struct AbsFMPTemplate
    w::Vector{Float64}                            # Primary Generator decision variables.               (GIVEN CONSTANT)
    gamma::Vector{Float64}                        # the bound for the demands on all the buses during a specific time  (GIVEN CONSTANT)
    q::Vector{Vector{Int}}                        # the secondary discrete decision variables           (GIVEN CONSTANT)
    d_hat::Vector{Float64}                        # the average testing demands vector.                 (GIVEN CONSTANT)

    v::Vector{Vector{VariableRef}}                # The slack decision variables for each of the previous demands.
    u::Vector{Vector{VariableRef}}                # The secondary continuous decision variables.
    # d::Vector{VariableRef}                        # The demand decision variable, as a giant vector.
    eta::VariableRef                              # The eta lower bound for all feasibility.
    lambda::Vector{Vector{VariableRef}}           # the dual decision variables.
    
    k::Int                                        # The iteration number from the ccga.

    sparse_vee::Bool
    demands_idx::Set{Int}                              # subset of indices indicating the demands constraints. 

end


"""
    Given the primary parameters and the secondary discrete decision variables
    meshed into a giant feasible set of many many constraints, we are
    searching for adversarial demands that can break our delivery system.
        * Determines the upper bound for the feasibility slack.

"""
@AbsFMPTemplate @ProblemTemplate mutable struct FMP <: AbsFMP
    
    rho_plus::Vector{VariableRef}                         # binary decision variables for the bilinear demands
    rho_minus::Vector{VariableRef}                        # binary decision variables for the bilinear demands
    xi_plus::Vector{Containers.DenseAxisArray}            # The binary decision variable for the extreme demands, vertices of the demands cube. 
    xi_minus::Vector{Containers.DenseAxisArray}           # The binary decision variable for the extreme demands. 
    
   
    function FMP()
        @warn("This constructor shouldn't be called except for debugging purposes. ")
        d̂ = 100.0*ones(size(MatrixConstruct.H, 2))
    return FMP(zeros(size(MatrixConstruct.B, 2)), 10.0, d̂) end

    function FMP(
        w::Vector,
        gamma::Vector,
        d_hat::Vector,
        model::Model=Model(HiGHS.Optimizer);
        sparse_vee::Bool=false
    )
        this = new()
        # Verify dimensions: 
        @assert length(w) == size(MatrixConstruct.B, 2) "w is having the wrong dimension when it's passed to the FMP. "
        @assert length(gamma) == MatrixConstruct.CONST_PROBLEM_PARAMETERS.HORIZON "gamma, the demand interval vector should be the same length as the time horizon, but it's not. " 
        this.w = w
        this.gamma = gamma
        this.q = Vector{Vector}()
        this.v = Vector{Vector}()
        this.u = Vector{Vector}()
        this.lambda = Vector{Vector}()
        this.k = 1
        this.model = model
        this.d_hat = d_hat
        this.con = Vector{Vector}()
        this.xi_plus = Vector{Containers.DenseAxisArray}()
        this.xi_minus = Vector{Containers.DenseAxisArray}()
        this.sparse_vee = sparse_vee

        IntroduceVariables!(this)
        PrepareConstraints!(this)
        PrepareObjective!(this)
    return this end

end


# The adding of the constraints API methods for FMP
"""
    THIS METHOD IS FOR POLYMORPHISM!!!
    It's being inherited by subtypes of AbsFMP. 
"""
function IntroduceVariables!(
    ::Type{AbsFMP}, 
    this::AbsFMP, 
    q_given::Union{Nothing, Vector{Float64}}
)
    k = this.k

    if this.sparse_vee
        starting, ending = MatrixConstruct.CON_GROUPS["Demand Balance"]
    else
        starting, ending = (1, size(MatrixConstruct.H, 1))
    end
    this.demands_idx = D = Set(starting:ending)
    D̃ = setdiff(Set(1:size(MatrixConstruct.H, 1)), D)|>collect|>sort

    # Prepare variable u for the model.
    push!(this.u, PrepareVariablesForTheModel!(this|>GetModel, :u, k))
    
    # Prepare for variable q, depending on whether a `q_given` is defined. 
    if q_given === nothing
        Decisionvariables = Vector()
        q = MatrixConstruct.q
        push!(Decisionvariables, zeros(Int, size(q[1])...)...)
        push!(Decisionvariables, zeros(Int, size(q[2])...)...)
        push!(Decisionvariables, zeros(Int, size(q[3])...)...)
        push!(this.q, Decisionvariables)
        this.eta = @variable(this.model, η >= 0)
    else
        push!(this.q, q_given)
    end

    # Prepare for variable v, the violations for all k. 
    v = @variable(
        this.model, 
        [1:length(MatrixConstruct.h)], 
        lower_bound=0, 
        base_name="v[$(k)]"
    )
    push!(this.v, v[:])
    
    if this.sparse_vee
        set_upper_bound.(v[D̃], 0)
    end

    # Parepare the variable lambda, the dual decision variables.
    λ = @variable(
        this.model, 
        [1:length(MatrixConstruct.h)], 
        upper_bound=0,lower_bound=-1,
        base_name="λ[$(k)]"
    )

    if this.sparse_vee
        delete_lower_bound.(λ[D̃])
    end
    
    push!(this.lambda, λ[:])

return this end



"""
    Prepare a new set of decision variables for the given instance of the problem.
        * u[k]::The continuous decision variable, for the kth iterations of the CCGA Algorithm.
            * u[k] will be a compositive decision variables for all the modeling decision variables in the original 
            problem.
        * v, λ follows a similar token.
    NOTE:
        * ALWAYS, run PrepareConstraints After a new constant "q[k]" is introduced to the system. 
        * This method is for super type AbsFMP and it's sharing it downwards. 
"""
function IntroduceVariables!(this::FMP, q_given::Union{Nothing, Vector{Float64}}=nothing)
    @assert q_given !== nothing || this.k == 1 "The conditions \"if q is nothing, then k == 1 for FMP is not true. \""
    IntroduceVariables!(AbsFMP, this, q_given)
    # Variables that we are using for the reformulations for the bi-linear constraints via the corner point assumptions

    model = this.model
    k = this.k
    B̄ = size(MatrixConstruct.H, 2)
    IdxTuples = findall(!=(0), MatrixConstruct.H).|>Tuple
    push!(this.xi_plus, @variable(model, [IdxTuples], base_name="Ξ⁺[$(k)]"))
    push!(this.xi_minus, @variable(model, [IdxTuples], base_name="Ξ⁻[$(k)]"))

    # Demand corner points decision variables: 
    if k == 1
        this.rho_plus = @variable(model, [1:B̄], Bin, base_name="ρ⁺")
        this.rho_minus = @variable(model, [1:B̄], Bin, base_name="ρ⁻")
    end

return this end




"""
    Prepare the constraints for the FMP problem, however, it will only add for all the most recent variables with
    r = k.
        * The class instance it keeping track of the interations count.
"""
function PrepareConstraints!(this::FMP)

    function addConstraints!(cons::JuMP.ConstraintRef)
        push!(this.con, cons)
    end

    k = this.k
    η = this.eta
    u = this.u[k]
    v = this.v[k]
    w = this.w
    H = MatrixConstruct.H
    B = MatrixConstruct.B
    G = MatrixConstruct.G
    C = MatrixConstruct.C
    h = MatrixConstruct.h
    λ = this.lambda[k]
    # d = this.d
    d̂ = this.d_hat
    γ = this.gamma
    #TODO: Change here for gamma
    Γ = kronecker(γ|>Diagonal, ones(MatrixConstruct.B̄)|>Diagonal)|>collect|>Diagonal 
    q = this.q[k]
    
    
    IdxTuples = findall(!=(0), H).|>Tuple
    model = this.model

    
    @constraint(model, η <= sum(v), base_name="eta objective: [$(k)]").|>addConstraints!
    # Dual constraints
    @constraint(model, λ'*C .<= 0, base_name="dual constraints: [$(k)]").|>addConstraints!

    # Sparse bilinear constraints setup with corner point assumptions:
    ρ⁺ = this.rho_plus
    ρ⁻ = this.rho_minus
    ξ⁺ = this.xi_plus[k]
    ξ⁻ = this.xi_minus[k]

    for (j, b) in IdxTuples
        @constraint(model, -ρ⁺[b] <= ξ⁺[(j, b)])|>addConstraints!
        @constraint(model, ξ⁺[(j, b)] <= ρ⁺[b])|>addConstraints!
        @constraint(model, λ[j] - (1 - ρ⁺[b])<=ξ⁺[(j, b)])|>addConstraints!
        @constraint(model, ξ⁺[(j, b)] <= λ[j] + (1 - ρ⁺[b]))|>addConstraints!
        @constraint(model, -ρ⁻[b] <= ξ⁻[(j, b)])|>addConstraints!
        @constraint(model, ξ⁻[(j, b)] <= ρ⁻[b])|>addConstraints!
        @constraint(model, λ[j] - (1 - ρ⁻[b]) <= ξ⁻[(j, b)])|>addConstraints!
        @constraint(model, ξ⁻[(j, b)] <= λ[j] + (1 - ρ⁻[b]))|>addConstraints!
    end
    @constraint(model, [b=1:size(H, 2)], ρ⁺[b] + ρ⁻[b] == 1, base_name="bin con:[$(k)]").|>addConstraints!

    # Bi-linear constraints
    
    @constraint(
        model,
        dot(λ, -H*d̂) + dot(-(H*Γ)[H.!=0], ξ⁺[:] .- ξ⁻[:]) + dot(λ, h - B*w - G*q) == sum(v),
        base_name="bilinear obj:[$(k)]"
    ).|>addConstraints!
    
    # new added demand feasibility constraints, the primal
    @constraint(
        model,
        B*w + G*q + C*u  + H*d̂ +  H*Γ*(ρ⁺ - ρ⁻) - v .<= h,
        base_name="opt con: [$(k)]"
    ).|>addConstraints!

    # Additional constraints for the sparse vee conditions that might get applied here.
    @constraint(model, sum(v) <= 500000)

return this end


"""
    Adversarial demands are coming from vertices of the hyper cube, this will recover
    the term (d̂ + γ.*ρ⁺ - γ.*ρ⁻), which represents the extreme demands which breaks the system.
        * It's used for the feasibility cut for the master problem to determine primary decision variable w̄.
        * Returns the decision variable from the model, explicit conversion to value vectors is needed.
"""
function GetDemandVertex(this::FMP)
    ρ⁺ = this.rho_plus.|>value 
    ρ⁻ = this.rho_minus.|>value
    γ = this.gamma.|>value
    #TODO: Change here for gamma
    Γ = kronecker(γ|>Diagonal, ones(MatrixConstruct.B̄)|>Diagonal)|>collect|>Diagonal
    d̂ = this.d_hat
return d̂ + Γ*(ρ⁺ - ρ⁻) end


"""
    Get the most recent ρ⁺, the worse demands for all the secondary configurations.
"""
function GetRhoPlus(this::FMP)
return this.rho_plus.|>value end


"""
    Get the most recent ρ⁻, the worst demands for all the secondary configurations.
"""
function GetRhoMinus(this::FMP)
return this.rho_minus.|>value end


"""
    Introduce a new secondary binary decision variable as a fixed variable for the system.
        * Introduces variable q[k] as a new configurations profile for the system. 
        * It will automatically assign the index, make the variables, and add the constraints for thsi specific secondary
        generator decision variables. 
"""
function Introduce!(this::FMP, q::Vector{Float64})
    this.k += 1
    @assert length(q) == size(MatrixConstruct.G, 2) "Wrong size for the passed in discrete secondary decision variables: q. "
    IntroduceVariables!(this, q)
    PrepareConstraints!(this)
return this end

"""
    Prepare the objective for the FMP problem, it's the sum of all v for all system configurations. 
"""
function PrepareObjective!(this::AbsFMP)
    @objective(this.model, Max, this.eta)
return this end

### ====================================================================================================================
### McCormickFMP: 
###     A convex relaxations for the collinear constraints appeared in the FMP strong duality objective. It uses the 
###     McCormick Envelope to achieve it. 
### ====================================================================================================================


@AbsFMPTemplate @ProblemTemplate mutable struct McCormickFMP <: AbsFMP
    
    demands_con_idx::Set{Int}                  # The set of indices that corresponds to the demand balance constraints. 
    not_demands_con_idx::Set{Int}              # The set of indices that is not the demand balance constraints. 
    d::Vector{VariableRef}
    xi::Vector{Containers.DenseAxisArray}

    function McCormickFMP(sparse_vee::Bool=true)
        @warn("This construction only exists for testing purposes!")
        return McCormickFMP(
            Model(HiGHS.Optimizer),
            40*ones(size(MatrixConstruct.H, 2)),
            zeros(size(MatrixConstruct.B, 2)),
            10; sparse_vee
        )
    end

    function McCormickFMP(
        model::Model, 
        d_hat::Vector, 
        w_bar::Vector, 
        gamma_bar; 
        sparse_vee::Bool=true
    )
        @error "The struct: McCormickFMP is DEPRECATED and shouldn't be used."
        this = new()
        this.w = w_bar
        this.gamma = gamma_bar
        this.q = Vector{Vector}()
        this.v = Vector{Vector}()
        this.u = Vector{Vector}()
        this.lambda = Vector{Vector}()
        this.k = 1
        this.model = model
        this.d_hat = d_hat
        this.con = Vector{Vector}()
        this.xi = Vector{Containers.DenseAxisArray}()
        this.sparse_vee = sparse_vee
        IntroduceVariables!(this)
        PrepareConstraints!(this)
        PrepareObjective!(this)
    return this end

end


"""
    Introduce all the decision variables for the McCormic version of the FMP. 
        * The upper bouned and lower bound for the variables are introduced as constraints in this function 
        as well. 
"""
function IntroduceVariables!(this::McCormickFMP, q_given::Union{Nothing, Vector{Float64}}=nothing)
    IntroduceVariables!(AbsFMP, this, q_given)
    γ̄ = this.gamma
    d̂ = this.d_hat
    k = this.k
    
    # The decision variable for McCormic Envelope. 
    # the ξ is Sparse! 
    IdxTuples = findall(!=(0), MatrixConstruct.H).|>Tuple
    push!(this.xi, @variable(this.model, [IdxTuples], base_name="ξ[$(k)]"))
    if k == 1
        d = this.d = @variable(this.model, d[1:size(MatrixConstruct.H, 2)])
        push!(
            this.con, @constraint(this.model, -γ̄ .<= d - d̂ .<= γ̄)...
        )
    end
return this end


"""
    Prepare the constraints for the FMP McCormic Formulations. 
"""
function PrepareConstraints!(this::McCormickFMP)
    function addConstraints!(cons::JuMP.ConstraintRef)
        push!(this.con, cons)
    end
    k = this.k
    η = this.eta
    u = this.u[k]
    v = this.v[k]
    w = this.w
    H = MatrixConstruct.H
    B = MatrixConstruct.B
    G = MatrixConstruct.G
    C = MatrixConstruct.C
    h = MatrixConstruct.h
    λ = this.lambda[k]
    ξ = this.xi[k]
    d = this.d
    d̂ = this.d_hat
    γ̄ = this.gamma
    q = this.q[k]
    # B̄ = size(H, 2)
    IdxTuples = findall(!=(0), H).|>Tuple
    model = this.model

    
    # Dual constraints
    @constraint(model, η <= sum(v), base_name="Objective lower Bound, k=$(k)").|>addConstraints!
    @constraint(model, λ'*C .<= 0, base_name="Dual λ'C: k=$(k)").|>addConstraints!
    
    # sparse Vee related constraints
    # @constraint(model, -1 .<= λ[D] .<= 0, base_name="λ range: k=$(k)").|>addConstraints!

    # Bi-linear equality constraints: The strong duality via the McCormic Envelope
    @constraint(
        model,
        dot(-H[H.!=0], ξ[:]) + dot(λ,h - B*w - G*q) == sum(v),
        base_name="strong duality: k=$(k)"
    ).|>addConstraints!
    for (j, b) in IdxTuples
        if j in this.demands_idx
            # (4a)
            @constraint(
                model,
                ξ[(j, b)] >= (d̂[b] - γ̄)*λ[j] - d[b] + d̂[b] - γ̄, 
                base_name="MC 4a D: k=$(k)"
            ).|>addConstraints!
            # (4c)
            @constraint(
                model, 
                ξ[(j, b)] <= (d̂[b] + γ̄)*λ[j] - d[b] + d̂[b] + γ̄, 
                base_name="MC 4c D: k=$(k)"
            ).|>addConstraints!
        else
            @constraint(
                model, 
                ξ[(j, b)] == 0, 
                base_bame="MC 4c ~D: k=$(k)"
            ).|>addConstraints!
        end
        # (4b)
        @constraint(
            model, 
            ξ[(j, b)] >= λ[j]*(d̂[b] + γ̄),
            base_name="MC 4b: k=$(k)"
        ).|>addConstraints!
        # (4d)
        @constraint(
            model, 
            ξ[(j, b)] <= λ[j]*(d̂[b] - γ̄),
            base_name="MC 4d: k=$(k)"
        ).|>addConstraints!
    end
    
    # operational constraints: 
    @constraint(
        model,
        B*w + G*q + C*u  + H*d - v .<= h,
        base_name="Operational constraints:k=$(k)"
    ).|>addConstraints!
    

return this end

"""
    Get the decision variable for demand as a vector of value. 
"""
function GetDemandVertex(this::McCormickFMP)
return this.d.|>value end


"""
    Add a new q variable for the FMP, and the q should be a solution from the FSP. 
"""
function Introduce!(this::McCormickFMP, q::Vector{Float64})
    this.k += 1
    @assert length(q) == size(MatrixConstruct.G, 2) "Wrong size for the passed in discrete secondary decision variables: q. "
    IntroduceVariables!(this, q)
    PrepareConstraints!(this)
return this end


### ====================================================================================================================
### GENETRIC FUNCTIONS FOR ALL ABOVE TYPES and Subsets of These types. 
###     - Generics Functions for abstract type of problems that contain a JuMP inner model, specific overrides too.
### ====================================================================================================================


"""
    Solve the inner JuMP optimization problem by directly invoking `optimize!` in JuMP on the model.
"""
function Solve!(this::Problem)
return optimize!(this|>GetModel) end


"""
    Get the JuMP instance from this model.
"""
function GetModel(this::Problem)
return this.model end



"""
    Produce a report for the MP (master problem), which also checks feasibility of the original problem and stuff.
        * Print out the string representation of the model.
        * Print out the constraints for the model.
        * Print out the objective for the model.
"""
function DebugReport(this::Problem, filename="debug_report")
    open("$(filename).txt", "w") do io
        write(io, this|>GetModel|>repr)
        write(io, "\n")
        this.con.|>(to_print) -> println(io, to_print)
        if !(this|>objective_value|>isnan)
            ["$(t[1]) = $(t[2])" for t in zip(this|>all_variables, this|>all_variables.|>value)].|>(to_print) -> println(io, to_print)
        end
    end
return this end


function PrintConstraintsGroup(con::Vector{JuMP.ConstraintRef}, filename::String="Constraints_report")
    open("$(filename).txt", "w") do io
        con.|>repr.|>(to_print) -> println(io, to_print)
    end
return end


"""
    Port out an variable for the user to write a functional to access a specific variable. 
        * Port out the decision variable from the model, if you give it the
            * Problem 
            * The symbol for the field for the problem. 
    User's Functional Should: 
        * Accept on parameters of type `JuMP.VariableRef`, and perform the wanted modificatons 
        for the decision variable. 
"""
function PortOutVariable!(port_out::Function, this::Problem, variable::Symbol)
    if variable in this|>typeof|>fieldnames
        return port_out(getfield(this, variable))
    end
    return (this |> GetModel)[variable]
return end


"""
    * Port out the model of the system as the first parameters of the context function. 
    * Port out the field of the instance as a dictionary of symbols ->  Field references as the second 
    parameters for the function. 
"""
function PortOutEverything!(port_out::Function, this::Problem)
    FieldDict = Dict([field => getfield(this, field) for field in fieldnames(this)])
return port_out(GetModel(this), FieldDict) end


# ----------------------------------------------------------------------------------------------------------------------

"""
    Get the value of the decision variable q, it's suitable for the type: Union{FSP, MP, FMP}. 
"""
function Getq(this::Union{FSP, MP})
return this.q.|>value.|>(x) -> round(x, digits=1)  end

"""
    Get the value of the u decision variables. 
"""
function Getu(this::Union{FSP, MP})
return this.u.|>value end


function Getv(this::Union{FSP, MP})
return this.v.|>value end


"""
    If the instance is an FMP, then it returns the variable v that is the most recent. 
"""
function Getv(this::FMP)
return this.v[end].|>value end



### ====================================================================================================================
### Methods Shared with JuMP.jl Model type
### ====================================================================================================================
 

"""
    * NaN is returned if there is no solution to the problem, which means it's either infeasible
    or unbounded.
"""
function objective_value(this::Problem)
    model = GetModel(this)
    status = primal_status(model)
    if status == NO_SOLUTION
        return NaN
    end
return JuMP.objective_value(model) end


"""
    List out the variable references for all the decision variables that are in the problem's model

"""
function all_variables(this::Problem) return this|>GetModel|>JuMP.all_variables end


"""
    Inherit the same indexer from the JuMP.Model. 
"""
function Base.getindex(this::Problem, param::Any)
return GetModel(this)[param] end


@info "CCGA Modeling tools has been loaded. "


