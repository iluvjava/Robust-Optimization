### ====================================================================================================================
# FSP: Lower bound searcher!
#   * Takes in a demands and it tries to satisfies it by choosing the secondary.
#   * discrete and contiuous variables.
#   * Templates are in the preliminaries. 
### ====================================================================================================================
"""
    Given demand from FMP, test how feasible it is to determine an Lower Bound
    for the feasibility slack variables.
"""
@ProblemTemplate mutable struct FSP <: Problem
    # model::Model
    u::Vector{VariableRef}
    v::VariableRef
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
        @warn("Sparse_vee option for FMP has been deprecated for the FSP, the option is now useless. ")
        this = new()
        this.model = model
        # this.gamma = gamma
        this.w = w
        this.d_star = d_star
        this.con = Vector{JuMP.ConstraintRef}()
        this|>IntroduceVariables!
        this|>AddConstraints!
        @objective(this.model, Min, this.v)
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
        this.v = @variable(this.model, v, lower_bound=0)

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
            @constraint(this.model, C*u + H*d + B*w + G*q .- v .<= h)...
        )
    return end


    function Base.getindex(this::FSP, index::Symbol)
    return this.model[index] end

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
        @assert length(gamma) == size(MatrixConstruct.H, 2) "gamma, the demand interval vector should be the same length as the time horizon, but it's not. " 
        @warn("Sparse Vee Feature has been DEPRECATED, this option is now useless. ")
        this.w = w
        this.gamma = gamma
        this.q = Vector{Vector}()
        this.v = Vector()   
        this.u = Vector{Vector}()
        this.lambda = Vector{Vector}()
        this.k = 1
        this.model = model
        this.d_hat = d_hat
        this.con = Vector{Vector}()
        this.xi_plus = Vector{Containers.DenseAxisArray}()
        this.xi_minus = Vector{Containers.DenseAxisArray}()
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

    # if this.sparse_vee
    #     starting, ending = MatrixConstruct.CON_GROUPS["Demand Balance"]
    # else
    #     starting, ending = (1, size(MatrixConstruct.H, 1))
    # end

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

    # TODO: CHANGE HERE
    v = @variable(
        this.model, 
        lower_bound=0, 
        base_name="v[$(k)]"
    )
    push!(this.v, v)
    
    

    # Parepare the variable lambda, the dual decision variables.
    # TODO: CHANGE HERE
    λ = @variable(
        this.model, 
        [1:length(MatrixConstruct.h)], 
        upper_bound=0, lower_bound=-1,
        base_name="λ[$(k)]"
    )

    
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
    Γ = γ|>Diagonal
    q = this.q[k]
    
    
    IdxTuples = findall(!=(0), H).|>Tuple
    model = this.model

    # TODO: CHANGE HERE
    @constraint(model, η <= v, base_name="eta objective: [$(k)]").|>addConstraints!
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
    # TODO: CHANGE HERE
    @constraint(
        model,
        dot(λ, -H*d̂) + dot(-(H*Γ)[H.!=0], ξ⁺[:] .- ξ⁻[:]) + dot(λ, h - B*w - G*q) == v,
        base_name="bilinear obj:[$(k)]"
    ).|>addConstraints!
    
    # new added demand feasibility constraints, the primal
    @constraint(
        model,
        (B*w + G*q + C*u  + H*d̂ +  H*Γ*(ρ⁺ - ρ⁻)) .- v .<= h,
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
    Γ = γ|>Diagonal
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
    @assert length(q) == size(MatrixConstruct.G, 2) "Wrong size for the passed in discrete secondary decision variables: q."
    
    IntroduceVariables!(this, q)
    PrepareConstraints!(this)
return this end

"""
    Prepare the objective for the FMP problem, it's the sum of all v for all system configurations. 
"""
function PrepareObjective!(this::AbsFMP)
    @objective(this.model, Max, this.eta)
return this end

