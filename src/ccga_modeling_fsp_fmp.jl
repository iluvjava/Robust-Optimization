### ====================================================================================================================
# FSP: Lower bound searcher!
#   * Takes in a demands and it tries to satisfies it by choosing the secondary.
#   * discrete and contiuous variables.
#   * Templates are in the preliminaries. 
### ====================================================================================================================
"""
Given demand vector from FMP that is within the uncertainty rainge, 
*FSP* tests how feasible it is, which determines an Lower Bound for the 
feasibility slack variable: `v`. 

### Fields

"""
@ProblemTemplate mutable struct FSP <: Problem
    # model::Model
    u::Vector{VariableRef}
    v::VariableRef
    q::Vector{VariableRef}

    # Given parameters
    w::Vector{Float64}
    d_star::Vector{Float64} 

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
    """
    function FSP(
        w::Vector{Float64}, 
        d_star::Vector{Float64}, 
        model::Model=Model(HiGHS.Optimizer)
    )   
        this = new()
        this.model = model
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
    lambda::Vector{Vector{VariableRef}}                   # The dual decision variable for each cut, lambda. 
    
   
    function FMP()
        @warn("This constructor shouldn't be called except for debugging purposes. ")
        d̂ = 100.0*ones(size(MatrixConstruct.H, 2))
    return FMP(zeros(size(MatrixConstruct.B, 2)), 10.0, d̂) end

    function FMP(
        w::Vector,
        gamma::Vector,
        d_hat::Vector,
        model::Model=Model(HiGHS.Optimizer);
    )
        this = new()
        # Verify dimensions: 
        @assert length(w) == size(MatrixConstruct.B, 2) "w is having the wrong dimension when it's passed to the FMP. "
        @assert length(gamma) == size(MatrixConstruct.H, 2) "gamma, the demand interval vector should be the same length as the time horizon, but it's not. " 
        this.w = w
        this.gamma = gamma
        this.k = 1
        this.model = model
        this.d_hat = d_hat
        this.q = Vector{Vector}()
        this.v = Vector()   
        this.u = Vector{Vector}()
        this.lambda = Vector{Vector}()
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
THIS METHOD IS FOR POLYMORPHISM!!!, It's being inherited by subtypes of AbsFMP. It sets up the following variables
for the abstract FMP types: 
* `q`: all zeros at the start, and then it's introduced by the FSP after the first iteration. If it's given, then it 
    will be assigned as a constant and stored in to the field of the type AbsFMP. 
* `u`: this is created as decision variable vector for each cut. 
    as a constant parameter in the field of the sub type. 
* `v`: The feasibility decision variable. 
* `η`: the lowerbound variable. 

### Argument: 
- `::Type{AbsFMP}`
- `this::AbsFMP`
- ` q_given::Union{Nothing, Vector{Float64}}`
"""
function IntroduceVariables!(
    ::Type{AbsFMP},
    this::AbsFMP, 
    q_given::Union{Nothing, Vector{Float64}}
)
    k = this.k

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
        lower_bound=0, 
        base_name="v[$(k)]"
    )
    push!(this.v, v)

return this end



"""
Prepare a new set of decision variables for the given instance of the *FMP* problem.
### NOTE:
* ALWAYS, run PrepareConstraints After a new constant "q[k]" is introduced to the system. 
* This method is for super type AbsFMP and it's sharing it downwards. 
### Arguments
- `this::FMP`
- `q_given::Union{Nothing, Vector{Float64}}=nothing`
"""
function IntroduceVariables!(this::FMP, q_given::Union{Nothing, Vector{Float64}}=nothing)
    @assert q_given !== nothing || this.k == 1 "The conditions \"if q is nothing, then k == 1 for FMP is not true. \""
    IntroduceVariables!(AbsFMP, this, q_given)
    # Variables that we are using for the reformulations for the bi-linear constraints via the corner point assumptions

    model = this.model
    k = this.k
    B̄ = size(MatrixConstruct.H, 2)
    # Parepare the variable lambda, the dual decision variables.
    λ = @variable(
        this.model, 
        [1:length(MatrixConstruct.h)], 
        upper_bound=0,
        base_name="λ[$(k)]"
    )
    IdxTuples = findall(!=(0), MatrixConstruct.H).|>Tuple
    push!(this.xi_plus, @variable(model, [IdxTuples], base_name="Ξ⁺[$(k)]"))
    push!(this.xi_minus, @variable(model, [IdxTuples], base_name="Ξ⁻[$(k)]"))
    
    push!(this.lambda, λ[:])
    # Demand corner points decision variables: 
    if k == 1
        this.rho_plus = @variable(model, [1:B̄], Bin, base_name="ρ⁺")
        this.rho_minus = @variable(model, [1:B̄], Bin, base_name="ρ⁻")
    end

return this end




"""
Prepare the constraints for the FMP problem, however, it will only add for all the most recent variables with
r = k.

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
    @constraint(
        model,
        dot(λ, -H*d̂) + dot(-(H*Γ)[H.!=0], ξ⁺[:] .- ξ⁻[:]) + dot(λ, h - B*w - G*q) == v,
        base_name="bilinear obj:[$(k)]"
    ).|>addConstraints!

    # Duality constraints
    @constraint(
        model, 
        sum(λ) >= -1,
        base_name="duality con:[$(k)]"
    ).|>addConstraints!

    # new added demand feasibility constraints, the primal
    @constraint(
        model,
        (B*w + G*q + C*u  + H*d̂ +  H*Γ*(ρ⁺ - ρ⁻)) .- v .<= h,
        base_name="opt con: [$(k)]"
    ).|>addConstraints!

    # Additional constraints for the sparse vee conditions that might get applied here.
    # @constraint(model, sum(v) <= 500000)

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
Introduce a cut figured out by the FSP instance. The cut is secondary decision variables. 
* Introduces variable q[k] as a new configurations profile for the system. 
* It will automatically assign the index, make the variables, and add the constraints for thsi specific secondary
generator decision variables. 
"""
function IntroduceCut!(this::FMP, q::Vector{Float64})
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


# ======================================================================================================================
# FMPH1, FMPH2: FMP with the bilinear heuristic. They are a "twin type". 
# ======================================================================================================================

"""
FMP with Bilinear search Heuristic. This is a super type of AbsFMP, it inherit all fields from it's super types 
recursively. 
### Fields
- `d::Vector{Float64}`: A fixed constant for the demands vector. 
### Constructor: 
- `w::Vector`
- `gamma::Vector`
- `d_hat::Vector`
"""
@AbsFMPTemplate @ProblemTemplate mutable struct FMPH1 <: AbsFMP
    
    d::Vector{Float64}   # decision variables: d. 
    lambda::Vector{Vector{VariableRef}}

    function FMPH1(
        w::Vector,
        gamma::Vector,
        d_hat::Vector,
        model::Model=Model(HiGHS.Optimizer)
    )
        this = new()
        this.w = w
        this.gamma = gamma
        this.k = 1
        this.model = model
        this.d_hat = d_hat
        this.q = Vector{Vector}()
        this.v = Vector()
        this.u = Vector{Vector}()
        this.lambda = Vector{Vector}()
        this.con = Vector{Vector}()
        # Construct the models 
        ort = [rand((-1, 1)) for __ in 1:size(MatrixConstruct.H, 2)]
        this.d = max.(ort.*this.gamma .+ this.d_hat, 0)
        IntroduceVariables!(this)
        PrepareConstraints!(this)
        PrepareObjective!(this)
        return this 
    end

end

"""
FMPH2 is the heuristic search for FMP where `lambda` is fixed to be a constant and the d vector is now a continuous 
    decision varaibles. 

### Fields
- `const_lambdas::Vector{Vector{Float64}}``
- `d::Vector{VariableRef}`

### Constructors
#### Args
- `w::Vector`
- `gamma::Vector`
- `d_hat::Vector`
- `lambdas::Vector{Vector}`
- `model::Model=Model(HiGHS.Optimizer)`
"""
@AbsFMPTemplate @ProblemTemplate mutable struct FMPH2 <: AbsFMP
    
    lambda::Vector{Vector{Float64}}
    d::Vector{VariableRef}

    function FMPH2(
        w::Vector,
        gamma::Vector,
        d_hat::Vector,
        lambda::Vector,
        model::Model=Model(HiGHS.Optimizer)
    )
        this = new()
        this.w = w
        this.gamma = gamma
        this.k = 1
        this.model = model
        this.d_hat = d_hat
        this.q = Vector{Vector}()
        this.v = Vector()
        this.u = Vector{Vector}()
        this.lambda = Vector()
        @assert length(lambda) == size(MatrixConstruct.H, 1)
        push!(this.lambda, lambda)
        this.con = Vector{Vector}()
        IntroduceVariables!(this)
        PrepareConstraints!(this)
        PrepareObjective!(this)
        return this
    end
end


"""
Introduce the JuMP variables to the FMPH model, no new variable are defined, all variables are the same as the abstract 
`IntroduceVariables!` method for the AbsFMP
### Arguments
- `this::Union{FMPH1, FMPH2}`: an instance of either FMP1, or FMP2. 
- `q_given::Union{Nothing, Vector{Float64}}=nothing`: An instance of `q_given`, either nothing whic is always 
that when k = 1 or something returned from the *FMP* when k is not 1. 
"""
function IntroduceVariables!(this::Union{FMPH1, FMPH2}, q_given::Union{Nothing, Vector{Float64}}=nothing)
    @assert q_given !== nothing || this.k == 1 "The conditions \"if q is nothing, then k == 1 for FMP is not true. \""
    IntroduceVariables!(AbsFMP, this, q_given)
    k = this.k
    if isa(this, FMPH1)
        # Parepare the variable lambda, the dual decision variables.
        λ = @variable(
            this.model, 
            [1:length(MatrixConstruct.h)], 
            upper_bound=0,
            base_name="λ[$(k)]"
        )
        push!(this.lambda, λ)
    else
        if k == 1
            this.d = @variable(this.model, d[1:size(MatrixConstruct.H, 2)], lower_bound=0)
        end
    end
    
    return
end

"""
Prepare the constraints for the FMPHDFixed type. Function should be used each time when a cut is introduce to the 
problem. 

"""
function PrepareConstraints!(this::Union{FMPH1, FMPH2})
    function AddConstraints!(cons::JuMP.ConstraintRef)
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
    d = this.d
    d̂ = this.d_hat
    γ = this.gamma
    Γ = γ|>Diagonal
    q = this.q[k]

    # eta lower bound constraint for each k
    @constraint(this.model, η <= v) .|> AddConstraints!
    # operatorational constraints for reach k
    @constraint(this.model, C*u .- v .<= H*d + h - B*w - G*q) .|> AddConstraints!
    # Duality constraints
    @constraint(this.model, λ'*C .<= 0) .|> AddConstraints!
    @constraint(this.model, sum(λ) >= -1) .|> AddConstraints!
    @constraint(this.model, dot(λ, H*d + h - B*w - G*q) == v) .|> AddConstraints!

    if isa(this, FMPH2) && k == 1
        @constraint(this.model, -this.gamma .<= this.d - this.d_hat .<= this.gamma) .|> AddConstraints!
    end

    return 
end

"""
Maximizes eta. 
"""
function PrepareObjective!(this::Union{FMPH1, FMPH2})
    model = this.model
    @objective(model, Max, this.eta)
    return this
end



"""
Introduce a new cut to the `FMPH1` with variables passed from the FSP. 
### Arguments
- `this::FMPHDFixed`: The type this function acts on. 
- `q::Vector{Float64}`: The q vector passed from the FSP instance. 
"""
function IntroduceCut!(this::FMPH1, q::Vector{Float64})
    @assert length(q) == size(MatrixConstruct.G, 2) "Wrong size for the passed in discrete secondary decision variables: q."
    this.k += 1
    IntroduceVariables!(this, q)
    PrepareConstraints!(this)
    return this 
end


"""
Introduce a new cut to the `FMPH2` with variables passed from the *FSP* and *FMPH1*
### Arguments
- `this::FMPHDFixed`: The type this function acts on. 
- `q::Vector{Float64}`: The q vector passed from the FSP instance. 
"""
function IntroduceCut!(this::FMPH2, lambda::Vector{Float64}, q::Vector{Float64})
    @assert size(MatrixConstruct.H, 1) == length(lambda)
    @assert length(q) == size(MatrixConstruct.G, 2) "Wrong size for the passed in discrete secondary decision variables: q."

    push!(this.lambda, lambda)
    this.k += 1
    IntroduceVariables!(this, q)
    PrepareConstraints!(this)
    return this 
end



"""
Try the binlinear search heuristic. 
"""
function TryHeuristic!(this::FMPH1, that::FMPH2, q::Vector{Float64})
    IntroduceCut!(this, q)
    Solve!(this)
    new_lambda = this.lambda[end].|>value  # new lambda from fmph1. 
    IntroduceCut!(that, new_lambda, q)
    Solve!(that)
    return that.eta|>value
end



function FirstHeuristic!(
    w::Vector, 
    gamma::Vector, 
    d_hat::Vector, 
    model::Model=Model(HiGHS.Optimizer)
)
    fmph1 = FMPH1(w, gamma, d_hat, model)
    Solve!(fmph1)
    fmph2 = FMPH2(w, gamma, d_hat, fmph1.lambda[1].|>value)
    Solve!(fmph2)
    
    return fmph2.eta|>value, fmph1, fmph2
end
