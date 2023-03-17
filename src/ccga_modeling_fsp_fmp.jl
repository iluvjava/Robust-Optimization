### ====================================================================================================================
# FSP: Lower bound searcher!
#   * Takes in a demands and it tries to satisfies it by choosing the secondary.
#   * discrete and contiuous variables.
#   * Templates are in the preliminaries. 
### ====================================================================================================================
"""
# Constructor
    FSP(w::Vector{Float64}, 
        d_star::Vector{Float64}, 
        model::Model=Model(HiGHS.Optimizer))

Given demand vector from FMP that is within the uncertainty rainge, *FSP* tests how feasible it is,
which determines an Lower Bound for the feasibility slack variable: `v`. 

# Fields

* constructed fields
    - `u::Vector{VariableRef}`: The continuous secondary decision variables. 
    - `v::VariableRef`: The slack variables for all the constraints. 
    - `q::Vector{VariableRef}`: The secondary discrete decision variable. 

* Given parameters
    - `w::Vector{Float64}`: The primary binary decision variables, from the master problem. 
    - `d_star::Vector{Float64} `: The adversarial demand suggested by the fmp algorithm. 
"""
@ProblemTemplate mutable struct FSP <: Problem

    "The continuous secondary decision variables. "
    u::Vector{VariableRef}
    "The slack variables for all the constraints, it's a scalar"
    v::VariableRef
    "The secondary discrete decision variable. "
    q::Vector{VariableRef}
    "The primary binary decision variables, from the master problem. "
    w::Vector{Float64}
    "The adversarial demand suggested by the fmp algorithm."
    d_star::Vector{Float64} 

    function FSP()
        @warn("This construction only exists for testing purposes!")
        return FSP(
            zeros(size(MatrixConstruct.B, 2)),
            10,
            ones(MatrixConstruct.d|>length)
        )
    end

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
    FMP(w::Vector,gamma::Vector,d_hat::Vector,model::Model=Model(HiGHS.Optimizer);)

Given the primary parameters and the secondary discrete decision variables
meshed into a giant feasible set of many many constraints, we are
searching for adversarial demands that can break our delivery system.

# Arguments 
- `w::Vector`: Where the length of the vector equals to the number of columns of B in the matrix construct module. 
This parameter should be determined by the master problem or the main problem, and it's a vector of binary value. 
- `gamma::Vector`: The uncertainty interval expressed as a vector. This gamma vector should be determined by the master problem. 
- `d_hat::Vector`: The center of the uncertainty interval. This vector should be none negative. 
- `model::Model=Model(HiGHS.Optimizer)`: The motimizer that FMP should be using. This model is an optimizer for solving the linear programming 
problem created for the FMP. 
"""
@AbsFMPTemplate @ProblemTemplate mutable struct FMP <: AbsFMP
    
    "binary decision variables for the bilinear demands: ρ⁺"
    rho_plus::Vector{VariableRef}
    "binary decision variables for the bilinear demands: ρ⁻"                        
    rho_minus::Vector{VariableRef}
    "The binary decision variable for the extreme demands, vertices of the demands cube: Ξ⁺"
    xi_plus::Vector{Containers.DenseAxisArray}
    "The binary decision variable for the extreme demands: Ξ⁻"
    xi_minus::Vector{Containers.DenseAxisArray}           
    "The dual decision variable for each cut, lambda: λ"
    lambda::Vector{Vector{VariableRef}}                   
    
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
    IntroduceVariables!(
        ::Type{AbsFMP},
        this::AbsFMP, 
        q_given::Union{Nothing, Vector{Float64}}
    )

THIS METHOD IS FOR POLYMORPHISM!!!, It's being inherited by subtypes of AbsFMP. It sets up the following variables
for the abstract FMP types: 
* `q`: all zeros at the start, and then it's introduced by the FSP after the first iteration. If it's given, then it 
    will be assigned as a constant and stored in to the field of the type AbsFMP. 
* `u`: this is created as decision variable vector for each cut. 
    as a constant parameter in the field of the sub type. 
* `v`: The feasibility decision variable. 
* `η`: the lowerbound variable. 

# Arguments
- `::Type{AbsFMP}`
- `this::AbsFMP`
- ` q_given::Union{Nothing, Vector{Float64}}`
"""
function IntroduceVariables!(
    ::Type{AbsFMP},
    this::AbsFMP, 
    q_given::Union{Nothing, Vector{N}}
) where {N <: Number}
    k = this.k

    # Prepare variable u for the model.
    push!(this.u, PrepareVariablesForTheModel!(this|>GetModel, :u, k))
    
    # Prepare for variable q, depending on whether a `q_given` is defined. 
    if q_given === nothing
        Decisionvariables = Vector()
        q = MatrixConstruct.q
        for q′ in q
            append!(Decisionvariables, zeros(Int, size(q′)...))
        end
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
    IntroduceVariables!(this::FMP, q_given::Union{Nothing, Vector{Float64}}=nothing)

Prepare a new set of decision variables for the given instance of the `FMP` problem.

# Arguments
- `this::FMP`
- `q_given::Union{Nothing, Vector{Float64}}=nothing`

# NOTE:
* ALWAYS, run PrepareConstraints After a new constant "q[k]" is introduced to the system. 
* This method is for super type AbsFMP and it's sharing it downwards. 

"""
function IntroduceVariables!(this::FMP, q_given::Union{Nothing, Vector{N}}=nothing) where {N<:Number}
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
Prepare the constraints for the initial FMP problem, without any cut. 
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
    IntroduceCut!(this::FMP, q::Vector{Float64}) 

Introduce a cut figured out by the FMP instance. The cut is secondary decision variables. 
---
* Introduces variable q[k] as a new configurations profile for the system. 
* It will automatically assign the index, make the variables, and add the constraints for thsi specific secondary
generator decision variables. 

# Arguments
- `this::FMP`: 
- `q::Vector{Float64}`: 
"""
function IntroduceCut!(this::FMP, q::Vector{Float64})
    this.k += 1
    @assert length(q) == size(MatrixConstruct.G, 2) "Wrong size for the passed in discrete secondary decision variables: q."
    IntroduceVariables!(this, q)
    PrepareConstraints!(this)
return this end

"""
    PrepareObjective!(this::AbsFMP)

Prepare the objective for the FMP problem, it's the sum of all v for all system configurations. 
"""
function PrepareObjective!(this::AbsFMP)
    @objective(this.model, Max, this.eta)
return this end


# ======================================================================================================================
# FMPH1, FMPH2: FMP with the bilinear heuristic. They are a "twin type". 
# ======================================================================================================================


"""
FMPH2 is the heuristic search for FMP where `lambda` is fixed to be a constant and the d vector is now a continuous 
    decision varaibles. 

# Fields
- `const_lambdas::Vector{Vector{Float64}}``
- `d::Vector{VariableRef}`

# Constructors Args
- `w::Vector`
- `gamma::Vector`
- `d_hat::Vector`
- `lambdas::Vector{Vector}`
- `model::Model=Model(HiGHS.Optimizer)`
"""
@ProblemTemplate @AbsFMPTemplate mutable struct FMPH2 <: AbsFMP
    
    "The dual for FMPH2, it's a constant. "
    lambda::Vector{Vector{Float64}}
    "The demand decision vector for FMPH2. "
    d::Vector{VariableRef}
    
    "The contain for the dualiy constraints where the constant lambds are involved here. "
    dual_cons::Vector{ConstraintRef}


    function FMPH2(
        w::Vector,
        gamma::Vector,
        d_hat::Vector,
        lambda::Vector,
        model::Model=Model(HiGHS.Optimizer)
    )
        @assert length(w) == size(MatrixConstruct.B, 2) "w is having the wrong dimension when it's passed to the FMPH2."
        @assert length(gamma) == size(MatrixConstruct.H, 2) "The uncertaity bound gamma vector is not having the compatible length with matrix B for FMPH2."
        @assert length(d_hat) == size(MatrixConstruct.H, 2) "The size of the d̂ is not compatible with the matrix for FMPH2. "
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
        # this.dual_con_idx = Vector{Int}()
        this.dual_cons = Vector{ConstraintRef}()
        IntroduceVariables!(this)
        PrepareConstraints!(this)
        PrepareObjective!(this)
        return this
    end
end


"""
FMP with Bilinear search Heuristic. This is a super type of AbsFMP, it inherit all fields from it's super types 
recursively. 
# Constructor
    FMPH1(
        w::Vector,
        gamma::Vector,
        d_hat::Vector,
        model::Model=Model(HiGHS.Optimizer);
        initial_demands::Union{Vector, Nothing}=nothing
    )

# Fields
- `d::Vector{Float64}`: A fixed constant for the demands vector. 
- `lambda::Vector{Vector{VariableRef}}`: The constants, lambda dual decision variable. It should be suggested 
by some instances of FMPH1. 
- `dual_con_idx::Vector{Int}`: A list of indices storing exactly where the constraints for the strong duality is inside 
of the constraint vector. 

# Constructor Positional Args: 
- `w::Vector`
- `gamma::Vector`
- `d_hat::Vector`
- `model::Model=Model(HiGHS.Optimizer)`

# Constructors Named Args: 
- `initial_demands::Union{Vector, Nothing}=nothing`: The initial demands with no cut added yet. 

"""
@ProblemTemplate @AbsFMPTemplate mutable struct FMPH1 <: AbsFMP
    
    d::Vector{Float64}   # decision variables: d. 
    lambda::Vector{Vector{VariableRef}}
   
    function FMPH1(
        w::Vector,
        gamma::Vector,
        d_hat::Vector,
        model::Model=Model(HiGHS.Optimizer);
        initial_demands::Union{Vector, Nothing}=nothing
    )
        @assert length(w) == size(MatrixConstruct.B, 2) "w is having the wrong dimension when it's passed to the FMPH1."
        @assert length(gamma) == size(MatrixConstruct.H, 2) "The uncertaity bound gamma vector is not having the compatible length with matrix B for FMPH1."
        @assert length(d_hat) == size(MatrixConstruct.H, 2) "The size of the d̂ is not compatible with the matrix for FMPH1. "
       
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
        if initial_demands === nothing
            ort = [1 for __ in 1:size(MatrixConstruct.H, 2)]
            this.d = max.(ort.*this.gamma .+ this.d_hat, 0)
        else
            @assert length(initial_demands) == size(MatrixConstruct.H, 2) "The size of the demands passed to FMPH2 is not in the right size."
            this.d = initial_demands
        end
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
function IntroduceVariables!(this::Union{FMPH1, FMPH2}, q_given::Union{Nothing, Vector{N}}=nothing) where {N<:Number}
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
    PrepareConstraints!(this::FMPH1)

Prepare the constraints for an instance of FMPH1. 
"""
function PrepareConstraints!(this::FMPH1)
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
    # Γ = γ|>Diagonal
    q = this.q[k]

    # eta lower bound constraint for each k
    @constraint(this.model, η <= v, base_name="eta con [$k]") .|> AddConstraints!
    # operatorational constraints for reach k
    @constraint(this.model, C*u + B*w + G*q + H*d .- v .<= h, base_name="opt con [$k]") .|> AddConstraints!
    # Duality constraints
    @constraint(this.model, λ'*C .<= 0, base_name="dual [$k]") .|> AddConstraints!
    @constraint(this.model, sum(λ) >= -1, base_name="dual [$k]") .|> AddConstraints!
    @constraint(this.model, dot(λ, -H*d + h - B*w - G*q) == v, base_name="strong dual [$k]") .|> AddConstraints!
   
    return this
end

"""
    PrepareConstraints!(this::FMPH2)

Prepare the constraints for an instance of FMPH2. 
"""
function PrepareConstraints!(this::FMPH2)
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
    
    d = this.d
    d̂ = this.d_hat
    γ = this.gamma
    # Γ = γ|>Diagonal
    q = this.q[k]

    # eta lower bound constraint for each k
    @constraint(this.model, η <= v, base_name="eta con [$k]") .|> AddConstraints!
    # operatorational constraints for reach k
    @constraint(this.model, C*u + B*w + G*q + H*d .- v .<= h, base_name="opt con [$k]") .|> AddConstraints!
    if length(this.lambda) >= k  
        λ = this.lambda[k]
        push!(this.dual_cons, @constraint(this.model, dot(λ, -H*d + h - B*w - G*q) == v, base_name="strong dual [$k]"))
    end
    if k == 1
        @constraint(this.model, -γ .<= this.d - d̂ .<= γ) .|> AddConstraints!
    end
    return this
end

"""
Get the demand decision variable by its value from the instance of FMPH2. If it's not yet solved, there will be an error. 
"""
function GetDemands(this::FMPH2)
    return this.d.|>value
end


"""
Get the value for all the decision variable lambds in the instance of FMPH1, if it's not yet solved, an error will occur. 
"""
function GetLambdas(this::FMPH1)
    return this.lambda.|>(l)->(l.|>value)
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
function IntroduceCut!(this::FMPH1, q::Vector{N}) where {N <: Number}
    @assert length(q) == size(MatrixConstruct.G, 2) "Wrong size for the passed in discrete secondary decision variables: q."
    this.k += 1
    IntroduceVariables!(this, q)
    PrepareConstraints!(this)
    return this 
end


"""
Introduce a new cut to the `FMPH2` with variables passed from the *FSP* and *FMPH1*
# Arguments
- `this::FMPHDFixed`: The type this function acts on. 
- `q::Vector{Float64}`: The q vector passed from the FSP instance. 
"""
function IntroduceCut!(this::FMPH2, q::Vector{N}) where {N <: Number}
    @assert length(q) == size(MatrixConstruct.G, 2) "Wrong size for the passed in discrete secondary decision variables: q."
    this.k += 1
    IntroduceVariables!(this, q)
    PrepareConstraints!(this)
    return this 
end


"""
    TryHeuristic!(this::FMPH1, that::FMPH2, q::Vector{Float64})

Try the binlinear search heuristic. Given the parameter for the FMP, return 2 instances of FMPH1, FMPH2, where 
FMPH1 has the demand variable as a constant randomly generated, and a q vector that is all zeros for both FMPH1, 2. 
FMPH1 solves for a value of lambda and then pass it to FMPH2 which then solves a demand vector.  

# Positional Arguments
- `this::FMPH1` 
- `that::FMPH2`
- `q::Vector{Float64}`
"""
function TryHeuristic!(this::FMPH1, that::FMPH2, q::Vector{Float64})
    IntroduceCut!(this, q)
    Solve!(this)
    IntroduceCut!(that, q)
    ChangeLambdas!(that, GetLambdas(this))
    Solve!(that)
    return this|>objective_value, that|>objective_value
end


"""
    TryHeuristic!(this::FMPH1, that::FMPH2)

Perform one exchange of the heurstic value. FMPH1 solves the system to obtain a lambda value, and pass the 
lambda value to FMPH2 to solve for the demand. 
"""
function TryHeuristic!(this::FMPH1, that::FMPH2)
    Solve!(this)
    ChangeLambdas!(that, GetLambdas(this))
    Solve!(that)
    return this|>objective_value, that|>objective_value
end

"""
    FirstHeuristic!(
        w::Vector, 
        gamma::Vector, 
        d_hat::Vector, 
        model1::Model=Model(HiGHS.Optimizer), 
        model2::Model=Model(HiGHS.Optimizer); 
        initial_demands::Union{Vector, Nothing}=nothing
    )

The FIRST Heuristic function accepts necessary parameters and then returns the instances of FMPH1, 2. We 
can also suggest good initial demands for it to speed things up. 

# Positional Arguments
- `w::Vector`
- `gamma::Vector`
- `d_hat::Vector`
- `model1::Model=Model(HiGHS.Optimizer)`: The empty model prepared for FMPH1. 
- `model2::Model=Model(HiGHS.Optimizer)`: The empty model prepared for FMPH2. 
# Keyword Arguments
- `initial_demands::Union{Vector, Nothing}=nothing`
"""
function FirstHeuristic!(
    w::Vector, 
    gamma::Vector, 
    d_hat::Vector, 
    model1::Model=Model(HiGHS.Optimizer), 
    model2::Model=Model(HiGHS.Optimizer); 
    initial_demands::Union{Vector, Nothing}=nothing
)
    fmph1 = FMPH1(w, gamma, d_hat, model1, initial_demands=initial_demands)
    Solve!(fmph1)
    fmph2 = FMPH2(w, gamma, d_hat, fmph1.lambda[1].|>value, model2)
    Solve!(fmph2)   
    return fmph1|>objective_value, fmph2|>objective_value, fmph1, fmph2
end


"""
    FirstHeuristic!(
        w::Vector, gamma::Vector, d_hat::Vector, 
        optimizer::DataType;
        kwargs...
    )

A variant of another function that only needs one function that produces a model instead of 2 separate models for 
FMPH1, 2. 
# Keyword Arguments
- `initial_demands::Union{Vector, Nothing}=nothing`: The initial demands to istantiate the instance of FMPH1. 
"""
function FirstHeuristic!(
    w::Vector, gamma::Vector, d_hat::Vector, 
    make_model_func::Function;
    kwargs...
)
return FirstHeuristic!(w, gamma, d_hat, make_model_func(), make_model_func(); kwargs...) 
end


"""
    ChangeLambdas!(this::FMPH2, lambdas::Vector{Vector{Float64}})

Given the lambds values for all cut, change the values of lambda in the underlying JuMP model using the new 
set of lambds, update the lambda references in this model. 
"""
function ChangeLambdas!(this::FMPH2, lambdas::Vector{Vector{Float64}})

    @assert this.k <= length(lambdas) "The set of lambdas suggested for FMPH2 is has length $(lambdas|>length), and the fmph2 has k=$(this.k). "
    H = MatrixConstruct.H    
    h = MatrixConstruct.h
    B = MatrixConstruct.B
    G = MatrixConstruct.G
    d = this.d
    w = this.w
    q = this.q
    v = this.v
    λ = lambdas
    
    # find the duality constraints and delete them. 
    for idx in 1:length(this.dual_cons)
        delete(this.model, this.dual_cons[idx])
    end
    this.dual_cons = Vector{ConstraintRef}()
    for k in 1:this.k
        push!(
            this.dual_cons, 
            @constraint(this.model, dot(λ[k], -H*d + h - B*w - G*q[k]) == v[k], base_name="strong dual [$k]")
        )
    end
    this.lambda = lambdas
    return this
end


"""
    BuildFMPH1(fmph1::FMPH1, d::Vector)

This function is going build another instance of the FMPH1 given the demands vector provided by FMPH2 
instance, and an existing instance of FMPH1 that we wish to replace the its demand vectors. 
------
# Arguments
- `fmph1::FMPH1`
- `d::Vector`
"""
function RebuildFMPH1(fmph1::FMPH1, d::Vector; model::Model=HiGHS.Optimizer)
    @assert length(d) == size(MatrixConstruct.H, 2) "The d used to build fmph1 is having a size of less than 2. "
    fmph1_new = FMPH1(fmph1.w, fmph1.gamma, fmph1.d_hat, model, initial_demands=d)
    for idx = 2:length(fmph1.q)
        IntroduceCut!(fmph1_new, fmph1.q[idx])
    end
    return fmph1_new
end


"""

# Constructor 

    FMPHStepper(w::Vector, gamma::Vector, d_hat::Vector, make_model_func::Function; kwargs...)

### key world arguments: 
    - `initial_demands::Union{Vector, Nothing}=nothing`: The initial demands for the instance of FMPH1. 

# Descriptions
A functor for executing the FMPH algorithm step by step and make collecting results simple and easy. 

"""
mutable struct FMPHStepper
    """
    An instance of the FMPH1, where d, demand is a constant and lambda, dual variable is a decision variable. 
    It will be updated throughout the algorithm 
    """
    fmph1::FMPH1
    """
    An instance of FMPH2, where d is the decision variable and lambdas are a constant, solved by FMPH1.
    It will be updated throughout the algorithm. 
    """
    fmph2::FMPH2
    "Objective value of fmph1 whenever is being solved. "
    obj1::Vector{Float64}
    "Objective value of the fmph2 whenever is being solved. "
    obj2::Vector{Float64}
    """
    A function that creates models for all of the FMPH1, FMPH2. This is a must because this struct 
    constantly build new instance of FMPH1, FMPH2. 
    """
    make_model_func::Function
    """
    An integer to keep track of the current execution cycle. One execution cycle refers to: 
    - Building up an instance of FMPH1 with the appropriate demand vector. 
    - Solving the instance of FMPH1 to get the dual variable λ
    - Putting the lambda into the instance of FMPH2 and solving that to get the demands. 
    """
    k::Int

    function FMPHStepper(w::Vector, gamma::Vector, d_hat::Vector, make_model_func::Function; kwargs...)
        this = new()
        obj1, obj2, fmph1, fmph2 = FirstHeuristic!(w, gamma, d_hat, make_model_func; kwargs...)
        this.fmph1 = fmph1
        this.fmph2 = fmph2
        this.obj1 = Vector{Float64}()
        this.obj2 = Vector{Float64}()
        push!(this.obj1, obj1)
        push!(this.obj2, obj2)
        this.make_model_func=make_model_func
        this.k = 0
        return this
    end

end


"""
    (this::FMPHStepper)()

Perform one step of the heuristic search and return the objective value of FMPH2. 
"""
function (this::FMPHStepper)()
    this.k += 1
    this.fmph1 = RebuildFMPH1(this.fmph1, GetDemands(this.fmph2); model=this.make_model_func())
    obj1, obj2 = TryHeuristic!(this.fmph1, this.fmph2)
    push!(this.obj1, obj1)
    push!(this.obj2, obj2)
    return this.fmph2|>objective_value
end


"""
    (this::FMPHStepper)(q::Vector{N}) where {N<:Number}

Introduce a cut to the FMPH1 and get the objective value of FMPH2 after the cut. This will perform 
one cycle of heuristic search. 
"""
function (this::FMPHStepper)(q::Vector{N}) where {N<:Number}
    this.k += 1
    this.fmph1 = RebuildFMPH1(this.fmph1, GetDemands(this.fmph2); model=this.make_model_func())
    obj1, obj2 = TryHeuristic!(this.fmph1, this.fmph2, q)
    push!(this.obj1, obj1) 
    push!(this.obj2, obj2)
    return this.fmph2|>objective_value
end


"""
Try a new random demands and rebuild the instance FMPH, and then perform one cycle of heuristic search using this new demands. 
The cuts will still be preserved. 
"""
function TryNewDemand(this::FMPHStepper, rand_strategy::Int=0) 
    gamma = this.fmph1.gamma
    d_hat = this.fmph1.d_hat
    if rand_strategy == 0   # Make a random demand vector. 
        ort = [1 for __ in 1:size(MatrixConstruct.H, 2)]  
        d = max.(ort.*gamma .+ d_hat, 0)
    else # mutate the current demands from fmph2 on 2 indices. Uniformally Randomly choose from the undertainty interval. 
        # d_old = fmph2|>GetDemands
        @error "Haven't been implemented yet. "
    end
    this.fmph1 = RebuildFMPH1(this.fmph1, d; model=this.make_model_func()) # rebuild the fmph1 using that new demand vector. 
    obj1, obj2 = TryHeuristic!(this.fmph1, this.fmph2)
    push!(this.obj1, obj1) 
    push!(this.obj2, obj2)
    return this.fmph2|>objective_value
end


"""
Get the demands vector's value from the current solved instance of FMPH2. 
"""
function GetDemands(this::FMPHStepper)
    return GetDemands(this.fmph2)
end


function objective_value(this::FMPHStepper)
    return this.fmph2|>objective_value
end
